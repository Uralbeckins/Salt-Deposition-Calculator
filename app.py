from typing import Dict
import os
import flet as ft
import pandas as pd
import math
import datetime

# Химические константы: молярные массы и заряды ионов
Mr = {'Cl': 35.453, 'SO4': 96.06, 'HCO3': 61.017, 'Ca': 40.078, 'Mg': 24.305,
      'NaK': 22.99, 'Ba': 137.327, 'Sr': 87.62}
z = {'Cl': -1, 'SO4': -2, 'HCO3': -1, 'Ca': 2, 'Mg': 2, 'NaK': 1, 'Ba': 2, 'Sr': 2}

# Параметры сульфатных солей: эмпирические коэффициенты и молярные массы
salts = {
    'Барит': {'a': 10.147, 'b': -4.946e-3, 'c': 11.650e-6, 'd': -5.315e-5,
              'e': -4.003, 'f': 2.787, 'g': -0.619, 'h': -1.850e-3,
              'ion': 'Ba', 'Mr': 233.387},

    'Целестин': {'a': 6.090, 'b': 2.237e-3, 'c': 5.739e-6, 'd': -4.197e-5,
                 'e': -2.082, 'f': 0.944, 'g': -8.650e-2, 'h': -1.873e-3,
                 'ion': 'Sr', 'Mr': 183.68},

    'Ангидрит': {'a': 2.884, 'b': 9.327e-3, 'c': 0.188e-6, 'd': -3.400e-5,
                 'e': -1.994, 'f': 1.267, 'g': -0.190, 'h': -3.195e-3,
                 'ion': 'Ca', 'Mr': 136.138},

    'Бассанит': {'a': 4.053, 'b': -1.792e-3, 'c': 11.400e-6, 'd': -7.070e-5,
                 'e': -1.734, 'f': 0.562, 'g': -2.170e-2, 'h': -6.436e-4,
                 'ion': 'Ca', 'Mr': 145.146},

    'Гипс': {'a': 3.599, 'b': -0.266e-3, 'c': 9.029e-6, 'd': -5.586e-5,
             'e': -0.847, 'f': 5.240e-2, 'g': 8.520e-2, 'h': -2.090e-3,
             'ion': 'Ca', 'Mr': 172.153}
}

# Вспомогательные преобразования единиц
def to_fahrenheit(T_c):
    return T_c * 9.0 / 5.0 + 32.0

def to_psi(P_atm):
    return P_atm * 14.6959

# Перевод мг/л в моль/л
def conc_molar(mg_per_l, Mr_value):
    # Обработка значений с запятой
    if isinstance(mg_per_l, str):
        mg_per_l = mg_per_l.replace(',', '.')
    mg_per_l = float(mg_per_l)
    return (mg_per_l / 1000.0) / Mr_value

# Расчет ионной силы раствора
def calc_I(c):
    return 0.5 * sum(c[i] * (z[i]**2) for i in c if i in z)


# Расчёт индекса насыщения карбонатной системы (кальцит)
def calc_SI_calcite(c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil):

    # Перевод температуры и давления
    TB, PB = to_fahrenheit(T), to_psi(P)
    I = calc_I(c)

    # Определение режима: есть газ, нет газа, либо информации недостаточно
    gas = False
    no_gas = False
    no_info = False

    if GF > 0 and Q > 0 and water_c > 0 and y_CO2 > 0:
        gas = True
    elif Q > 0 and water_c > 0 and X_CO2 > 0:
        no_gas = True
    else:
        no_info = True

    # Расчёт при наличии газа
    if gas == True:
        Qj, wc, GF, yCO2_prime = float(Q), float(water_c), float(GF), float(y_CO2)
        BWPD = Qj * wc / 100 * 6.28981
        BOPD = Qj * (100 - wc) / 100 * 6.28981
        MMCFD = Qj * (100 - wc) / 100 * GF * 35.3147 / 1e6

        tb05, pb05 = TB**0.5, PB**0.5
        fCO2 = math.exp(-7.66e-4 + 8.0e-6 * tb05 - 2.11e-9 * tb05 * pb05)

        denom = 1 + PB * fCO2 * (5*BWPD + 10*BOPD) * 1e-9 / ((TB + 460) * max(MMCFD, 1e-12))
        yCO2 = yCO2_prime / denom if denom != 0 else yCO2_prime
        denom2 = PB * yCO2 * fCO2 if PB * yCO2 * fCO2 != 0 else 1e-12

        return math.log10(max(c.get('Ca', 0) * c.get('HCO3', 0), 1e-30) / denom2) + \
               6.039 + 0.014463 * TB - 9.44e-5 * TB**2 - 6.185e-5 * PB - \
               1.895 * I**0.5 + 0.662 * I + 0.03654 * I**1.5 - 8e-5 * TB * I**0.5

    # Расчёт при отсутствии газа, но с данными по нефти
    if no_gas == True:
        Qj, wc, XCO2, Mroil, rhooil = float(Q), float(water_c), float(X_CO2), float(Mroil), float(dens_oil)
        Vw = Qj * wc / 100 * 1000
        Vo = Qj * (100 - wc) / 100 * 1000
        nCO2 = Vo * rhooil * 1000 / Mroil * XCO2 if Mroil != 0 else 0
        cCO2 = 7289 * nCO2 / (Vw + 3.04 * Vo) if (Vw + 3.04 * Vo) != 0 else 1e-12

        return math.log10(max(c.get('Ca', 0) * c.get('HCO3', 0), 1e-30) / cCO2) + \
               3.801 + 0.008115 * T - 9.028e-6 * T**2 - 7.419e-5 * P - \
               1.961 * I**0.5 + 0.695 * I - 0.01136 * I**1.5 - 1.604e-5 * T * I**0.5

    # Расчёт по упрощённой формуле при недостатке информации
    if no_info == True:
        pH = float(pH)
        return math.log10(max(c.get('Ca', 0) * c.get('HCO3', 0), 1e-30)) + \
               pH - 2.53 + 0.008943 * TB + 1.886e-6 * TB**2 - 4.855e-5 * PB - \
               1.470 * I**0.5 + 0.316 * I + 0.0537 * I**1.5 + 0.001297 * TB * I**0.5


# Индекс насыщения для сульфатных солей
def calc_SI_sulfate(salt, c, T, P):

    TB, PB = to_fahrenheit(T), to_psi(P)
    I = calc_I(c)

    # Равновесная константа взаимодействия Ме²⁺–SO₄²⁻
    logK_st = 2.301 + 0.00174 * TB + 4.553e-6 * TB**2 - 7.801e-6 * PB - \
              3.969 * I**0.5 + 2.280 * I - 0.459 * I**1.5 - 0.0006037 * TB * I**0.5
    K_st = 10**logK_st

    # Концентрации для расчёта свободного SO₄²⁻
    sum_cM = c.get('Ca', 0) + c.get('Mg', 0) + c.get('Ba', 0) + c.get('Sr', 0)
    temp = K_st * (sum_cM - c.get('SO4', 0))
    a = 1 + temp

    sqrt_term = math.sqrt(max(a**2 + 4*K_st*c.get('SO4', 0), 0))
    cSO4_free = (-a + sqrt_term) / (2*K_st) if K_st != 0 else c.get('SO4', 0)

    # Свободная концентрация соответствующего катиона
    ion_name = salts[salt]['ion']
    cM_free = c.get(ion_name, 0) / (1 + K_st*cSO4_free) if (1 + K_st*cSO4_free) != 0 else c.get(ion_name, 0)

    # Расчет равновесной константы конкретной соли
    coef = salts[salt]
    logK_c = coef['a'] + coef['b'] * TB + coef['c'] * TB**2 + coef['d'] * PB + \
             coef['e'] * I**0.5 + coef['f'] * I + coef['g'] * I**1.5 + coef['h'] * TB * I**0.5

    return math.log10(max(cM_free * cSO4_free, 1e-30)) + logK_c


# Максимально возможное осаждение конкретной соли (итеративный поиск)
def get_max_dc(salt, c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil):

    # Определяем пары ионов для соли
    cation, anion = ('Ca', 'HCO3') if salt == 'Кальцит' else (salts[salt]['ion'], 'SO4')

    # Максимально возможный расход ионов
    limit = max(0, min(c.get(cation, 0), c.get(anion, 0)))
    low, high = 0.0, limit

    # Начальный SI
    if salt == 'Кальцит':
        SI = calc_SI_calcite(c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil)
    else:
        SI = calc_SI_sulfate(salt, c, T, P)

    if SI is None or SI <= 0:
        return 0.0

    i = 0
    mid = 0

    # Двоичный поиск количества осадка
    while SI is not None and abs(SI) > 1e-6:

        mid = (low + high) / 2
        new_c = c.copy()

        new_c[cation] = max(0, new_c.get(cation, 0) - mid)
        new_c[anion] = max(0, new_c.get(anion, 0) - mid)

        if salt == 'Кальцит':
            SI = calc_SI_calcite(new_c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil)
        else:
            SI = calc_SI_sulfate(salt, new_c, T, P)

        if SI is not None and SI > 0:
            low = mid
        else:
            high = mid

        i += 1

    return mid


# Главный расчёт масс осадков с последовательным осаждением
def calculate_masses(c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil, N_steps=1000):

    # Стартовые значения SI
    SI_results = {'Кальцит': calc_SI_calcite(c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil)}

    for salt in salts:
        SI_results[salt] = calc_SI_sulfate(salt, c, T, P) if c.get('SO4', 0) > 0 else -100.0

    # Список пересыщенных солей
    active = [k for k, v in SI_results.items() if v is not None and v > 0]

    # Предварительное вычисление максимальных осаждений
    delta_max = {s: get_max_dc(s, c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil) for s in active}
    delta = {s: delta_max[s] / N_steps for s in active}

    precip = {s: 0.0 for s in active}
    current_c = c.copy()

    # Основной итерационный цикл
    for _ in range(N_steps):

        I_step = calc_I(current_c)
        to_precipitate = []

        # Выбираем соли, которые всё ещё пересыщены
        for s in active:
            si = calc_SI_calcite(current_c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil) if s == 'Кальцит' \
                 else calc_SI_sulfate(s, current_c, T, P)

            if si is not None and si > 0:
                to_precipitate.append(s)

        # Реальное осаждение с проверкой доступных ионов
        for s in to_precipitate:
            cation = 'Ca' if s == 'Кальцит' else salts[s]['ion']
            anion = 'HCO3' if s == 'Кальцит' else 'SO4'

            d_actual = min(delta[s], current_c.get(cation, 0), current_c.get(anion, 0))
            if d_actual <= 0:
                continue

            current_c[cation] -= d_actual
            current_c[anion] -= d_actual
            precip[s] += d_actual

    # Перевод осаждённых количеств в г/л
    masses = {}
    si_res = {}
    sum_masses = 0

    for s in precip:
        mr = 100.087 if s == 'Кальцит' else salts[s]['Mr']
        masses[s] = round(precip[s] * mr, 3)
        si_res[s] = format(precip[s], '.3e')
        sum_masses += masses[s]
    # Добавляем нули для солей, которых нет в осадках
    for key in ['Кальцит', *salts.keys()]:
        si_res.setdefault(key, 0.0)
        masses.setdefault(key, 0.0)

    masses['Общая масса солей'] = round(sum_masses, 4)
    si_res['Общая масса солей'] = '-'
    
    return masses, si_res


DOWNLOADS_DIR = os.path.join(os.path.expanduser("~"), "Downloads")


def main(page: ft.Page):
    # =============================================================================
    # НАСТРОЙКА СТРАНИЦЫ
    # =============================================================================
    table = None
    page.title = "Калькулятор осаждения солей"
    page.padding = 20
    page.window.width = 900
    page.window.height = 1200
    page.window.top = 0 
    page.scroll = ft.ScrollMode.AUTO
    page.window.maximized = True

    # =============================================================================
    # ФАЙЛОВЫЙ ВЫБОР
    # =============================================================================
    file_picker = ft.FilePicker()
    page.overlay.append(file_picker)

    # =============================================================================
    # ВВОД ПАРАМЕТРОВ: ТЕМПЕРАТУРА, ДАВЛЕНИЕ, pH
    # =============================================================================
    temp_input = ft.TextField(label="Температура, °C", label_style=ft.TextStyle(size=18), value="60", width=200)
    press_input = ft.TextField(label="Давление, атм", label_style=ft.TextStyle(size=18), value="30", width=200)
    ph_input = ft.TextField(label="pH", label_style=ft.TextStyle(size=18), value="7.42", width=200)
    
    # =============================================================================
    # ПАРАМЕТРЫ СКВАЖИНЫ
    # =============================================================================
    Q_input = ft.TextField(label="Дебит жидкости, м³/сут", label_style=ft.TextStyle(size=18), value="0", width=200)
    water_input = ft.TextField(label="Обводненность, %", label_style=ft.TextStyle(size=18), value="0", width=200)
    gas_input = ft.TextField(label="Газовый фактор, м³/м³", label_style=ft.TextStyle(size=18), value="0", width=200)

    # =============================================================================
    # ИОННЫЙ СОСТАВ ВОДЫ
    # =============================================================================
    cl_input = ft.TextField(label="Cl⁻, мг/л", label_style=ft.TextStyle(size=18), value="11838.8", width=200)
    so4_input = ft.TextField(label="SO₄²⁻, мг/л", label_style=ft.TextStyle(size=18), value="122", width=200)
    hco3_input = ft.TextField(label="HCO₃⁻, мг/л", label_style=ft.TextStyle(size=18), value="164.7", width=200)
    ca_input = ft.TextField(label="Ca²⁺, мг/л", label_style=ft.TextStyle(size=18), value="641.28", width=200)
    mg_input = ft.TextField(label="Mg²⁺, мг/л", label_style=ft.TextStyle(size=18), value="60.75", width=200)
    nak_input = ft.TextField(label="Na⁺ и K⁺, мг/л", label_style=ft.TextStyle(size=18), value="6892.18", width=200)
    ba_input = ft.TextField(label="Ba²⁺, мг/л", label_style=ft.TextStyle(size=18), value="1000", width=200)
    sr_input = ft.TextField(label="Sr²⁺, мг/л", label_style=ft.TextStyle(size=18), value="102", width=200)

    # =============================================================================
    # НЕФТЬ И CO2
    # =============================================================================
    oil_density_input = ft.TextField(label="Плотность нефти, кг/м³", label_style=ft.TextStyle(size=18), value="0", width=200)
    oil_mr_input = ft.TextField(label="Средняя молекулярная масса пластовой нефти, г/моль", label_style=ft.TextStyle(size=18), value="0", width=450)
    co2_in_gas_input = ft.TextField(label="Мольная доля CO₂ в газовой фазе", label_style=ft.TextStyle(size=18), value="0", width=300)
    co2_in_oil_input = ft.TextField(label="Мольная доля CO₂ в пластовой нефти", label_style=ft.TextStyle(size=18), value="0", width=300)
    
    # =============================================================================
    # КОЛОНКА ДЛЯ ВЫВОДА РЕЗУЛЬТАТОВ
    # =============================================================================
    results_column = ft.Column([], spacing=10, horizontal_alignment=ft.CrossAxisAlignment.START)

    # =============================================================================
    # ОБЯЗАТЕЛЬНЫЕ ПОЛЯ
    # =============================================================================
    required_fields = [temp_input, press_input, ph_input, ca_input, cl_input, so4_input, hco3_input, mg_input, nak_input, ba_input, sr_input]

    def on_field_change(e):
        field = e.control
        if field.value and field.value.strip():
            field.border_color = None
        else:
            field.border_color = ft.Colors.RED_400
        page.update()

    for f in required_fields:
        f.on_change = on_field_change

    # =============================================================================
    # НОРМАЛИЗАЦИЯ ЧИСЛОВЫХ ПОЛЕЙ (ЗАМЕНА ',' НА '.')
    # =============================================================================
    numeric_fields = [
        temp_input, press_input, ph_input,
        Q_input, water_input, gas_input,
        cl_input, so4_input, hco3_input, ca_input, mg_input, nak_input, ba_input, sr_input,
        oil_density_input, oil_mr_input, co2_in_gas_input, co2_in_oil_input
    ]

    def normalize_input(e):
        field = e.control
        if ',' in field.value:
            field.value = field.value.replace(',', '.')
            page.update()

    for f in numeric_fields:
        f.on_change = normalize_input

    # =============================================================================
    # ФУНКЦИЯ РАСЧЕТА ПРИ НАЖАТИИ КНОПКИ
    # =============================================================================
    def calculate_click(e):
        try:
            # Проверка обязательных полей
            required_fields = [temp_input, press_input]
            missing = [f.label for f in required_fields if not (f.value and f.value.strip())]
            if missing:
                for f in required_fields:
                    f.border_color = ft.Colors.RED_400 if not (f.value and f.value.strip()) else None
                page.open(ft.SnackBar(ft.Text("Заполните обязательные поля!"), bgcolor=ft.Colors.RED))
                page.update()
                return

            nonlocal table

            # --- Считываем параметры ---
            T = float(temp_input.value or 0)
            P = float(press_input.value or 0)
            pH = float(ph_input.value or 0)
            GF = float(gas_input.value or 0)
            Q = float(Q_input.value or 0)
            water_c = float(water_input.value or 0)
            y_CO2 = float(co2_in_gas_input.value or 0)
            X_CO2 = float(co2_in_oil_input.value or 0)
            Mroil = float(oil_mr_input.value or 0)
            dens_oil = float(oil_density_input.value or 0)
            
            # --- Ионы ---
            ions = {}
            if cl_input.value: ions['Cl'] = conc_molar(float(cl_input.value), Mr['Cl'])
            if so4_input.value: ions['SO4'] = conc_molar(float(so4_input.value), Mr['SO4'])
            if hco3_input.value: ions['HCO3'] = conc_molar(float(hco3_input.value), Mr['HCO3'])
            if ca_input.value: ions['Ca'] = conc_molar(float(ca_input.value), Mr['Ca'])
            if mg_input.value: ions['Mg'] = conc_molar(float(mg_input.value), Mr['Mg'])
            if nak_input.value: ions['NaK'] = conc_molar(float(nak_input.value), Mr['NaK'])
            if ba_input.value: ions['Ba'] = conc_molar(float(ba_input.value), Mr['Ba'])
            if sr_input.value: ions['Sr'] = conc_molar(float(sr_input.value), Mr['Sr'])
            c = ions

            # --- Расчет масс осадков ---
            masses, si_res = calculate_masses(c, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil)

            # --- Отображение результатов ---
            results_column.controls.clear()

            # Форматирование вывода: показывать 0 вместо 0.0
            def format_display(v):
                try:
                    if isinstance(v, (int, float)) and float(v) == 0.0:
                        return '0'
                    return str(v)
                except Exception:
                    return str(v)

            # Заголовок
            results_column.controls.append(
                ft.Container(
                    content=ft.Column([
                        ft.Text("Результаты расчета", size=22, weight=ft.FontWeight.BOLD),
                        ft.Text("Ожидаемые массы осадков:", size=16, weight=ft.FontWeight.BOLD)
                    ]),
                    padding=10,
                    alignment=ft.alignment.top_left
                )
            )

            # Таблица с результатами
            if table is None:
                columns = [
                    ft.DataColumn(ft.Text('Вещество')),
                    ft.DataColumn(ft.Text('Масса, г/л')),
                    ft.DataColumn(ft.Text('Индекс насыщения'))
                ]
                rows = [ft.DataRow([ft.DataCell(ft.Text(salt)), ft.DataCell(ft.Text(format_display(masses[salt]))), ft.DataCell(ft.Text(format_display(si_res[salt])))])
                    for salt in ['Кальцит', 'Барит', 'Целестин', 'Ангидрит', 'Бассанит', 'Гипс', 'Общая масса солей']]
                table = ft.DataTable(
                    columns=columns,
                    rows=rows,
                    border=ft.border.all(1),
                    border_radius=4,
                    horizontal_lines=ft.border.BorderSide(1, ft.Colors.BLACK),
                    vertical_lines=ft.border.BorderSide(1, ft.Colors.BLACK),
                    column_spacing=50,
                    sort_ascending=True
                )
            else:
                table.rows = [ft.DataRow([ft.DataCell(ft.Text(salt)), ft.DataCell(ft.Text(format_display(masses[salt]))), ft.DataCell(ft.Text(format_display(si_res[salt])))])
                              for salt in ['Кальцит', 'Барит', 'Целестин', 'Ангидрит', 'Бассанит', 'Гипс', 'Общая масса солей']]

            table_container = ft.Container(content=table, padding=10, alignment=ft.alignment.top_left)
            results_column.controls.append(table_container)
            page.update()

        except Exception as ex:
            sb = ft.SnackBar(content=ft.Text(f"Ошибка: {str(ex)}"), bgcolor=ft.Colors.RED)
            page.open(sb)
            sb.open = True
            page.update()

    # =============================================================================
    # СМЕНА ТЕМЫ
    # =============================================================================
    def toggle_theme(e):
        global nav_color 
        if page.theme_mode == ft.ThemeMode.LIGHT:
            page.theme_mode = ft.ThemeMode.DARK
            theme_btn.icon = ft.Icons.LIGHT_MODE
            nav_bar.bgcolor = ft.Colors.BLUE_GREY_900
            nav_bar.inactive_color = ft.Colors.GREY_500
            nav_bar.active_color=ft.Colors.CYAN_ACCENT_400
        else:
            page.theme_mode = ft.ThemeMode.LIGHT
            theme_btn.icon = ft.Icons.DARK_MODE
            nav_bar.bgcolor = ft.Colors.LIGHT_BLUE_300
            nav_bar.inactive_color = ft.Colors.GREY_700
            nav_bar.active_color = ft.Colors.BLUE_800

        page.update()

    # Кнопка для смены темы
    theme_btn = ft.IconButton(
                    icon = ft.Icons.DARK_MODE,
                    on_click = toggle_theme,
                    tooltip='Сменить тему'
                )

    # =============================================================================
    # ФУНКЦИИ ДЛЯ СТРАНИЦЫ РУЧНОГО ВВОДА
    # =============================================================================
    def home_page():
        return ft.Container(
            content=ft.Column([
                ft.Row(
                    [ft.Text("Калькулятор расчета рисков солеотложения по Оддо-Томсону", size=24, weight=ft.FontWeight.BOLD, text_align=ft.TextAlign.CENTER),
                    theme_btn],
                    alignment=ft.MainAxisAlignment.SPACE_BETWEEN 
                ),
                ft.Divider(height=20),

                # Термобарические условия
                ft.Container(
                    content=ft.Column([
                        ft.Text("Термобарические условия", size=18, weight=ft.FontWeight.BOLD),
                        ft.Row([temp_input, press_input], wrap=True)
                    ], spacing=20),
                    padding=10
                ),
                ft.Divider(),

                # Параметры скважины
                ft.Container(
                    content=ft.Column([
                        ft.Text("Показатели работы скважины", size=18, weight=ft.FontWeight.BOLD),
                        ft.Row([Q_input, water_input, gas_input], wrap=True)
                    ], spacing=20),
                    padding=10
                ),
                ft.Divider(),

                # Ионный состав
                ft.Container(
                    content=ft.Column([
                        ft.Text("Ионный состав воды", size=18, weight=ft.FontWeight.BOLD),
                        ft.Row([
                            ft.Column([ph_input, ca_input, ba_input]),
                            ft.Column([nak_input, mg_input, sr_input]),
                            ft.Column([hco3_input, so4_input, cl_input])
                        ], spacing=20, wrap=True)
                    ], spacing=10),
                    padding=10
                ),
                ft.Divider(),

                # Нефть и CO2
                ft.Container(
                    content=ft.Column([
                        ft.Text("Информация о нефти и CO2", size=18, weight=ft.FontWeight.BOLD),
                        ft.Row([oil_density_input], wrap=True),
                        ft.Row([oil_mr_input], wrap=True),
                        ft.Row([co2_in_gas_input, co2_in_oil_input], wrap=True)
                    ], spacing=20),
                    padding=10
                ),

                # Кнопка расчета
                ft.Container(
                    content=ft.ElevatedButton(
                        content=ft.Text("Рассчитать", size=20, weight=ft.FontWeight.NORMAL),
                        on_click=calculate_click,
                        width=400,
                        height=50,
                        style=ft.ButtonStyle(color=ft.Colors.WHITE, bgcolor=ft.Colors.BLUE)
                    ),
                    padding=10,
                    alignment=ft.alignment.center
                ),

                # Результаты
                results_column
            ]),
            padding=20,
        )

    # =============================================================================
    # ФУНКЦИЯ ДЛЯ СОЗДАНИЯ ФАЙЛА ШАБЛОНА
    # =============================================================================
    def create_template(e):
        columns = [
            "№", "Скважина", "Дата отбора",
            "Cl⁻, мг/л", "SO₄²⁻ мг/л", "HCO₃, мг/л", "Ca²⁺, мг/л", "Mg²⁺, мг/л", "Na⁺ и K⁺, мг/л", "Ba²⁺, мг/л", "Sr²⁺, мг/л",
            "рН", "Плотность воды, г/см³", "Температура, °С", "Давление, атм", "Плотность нефти, кг/м³",
            "Средняя молекулярная масса пл. нефти, г/моль", "Мольная доля CO₂ в газовой фазе", "Мольная доля CO₂ в пл. нефти",
            "Дебит жидкости, м³/сут", "Обводненность, %", "Газовый фактор, м³/м³"
        ]

        temp_df = pd.DataFrame(columns=columns)
        template_path = os.path.join(DOWNLOADS_DIR, "Шаблон.xlsx")
        temp_df.to_excel(template_path, index=False)
        page.open(ft.SnackBar(ft.Text(f"* Шаблон создан, проверьте папку \"Загрузки\""), bgcolor=ft.Colors.GREEN, duration=9000))
        page.update()

    calculate_excel_btn = ft.ElevatedButton(
        content=ft.Row([
            ft.Text("3. Рассчитать", size=16),
            ft.Icon(ft.Icons.CALCULATE)
        ], alignment=ft.MainAxisAlignment.CENTER, spacing=5),
        on_click=None,  # пока неактивна
        width=280,
        height=50,
        disabled=True,
        style=ft.ButtonStyle(
            color=ft.Colors.WHITE,
            bgcolor=ft.Colors.GREY  # неактивная серый
        )
    )

    
    preview_container=ft.Container(
                        content=ft.Text("Здесь появится таблица c результатами расчета", size=20, color=ft.Colors.GREY),
                        width=1200,
                        height=400,
                        alignment=ft.alignment.center,
                        border=ft.border.all(1, ft.Colors.GREY),
                        border_radius=30,
                        margin=ft.margin.only(top=10, bottom=20),
                        expand=False
                                )
    result_df = None

    # =============================================================================
    # ФУНКЦИЯ ДЛЯ ОБРАБОТКИ ФАЙЛА EXCEL И СОХРАНЕНИЯ РЕЗУЛЬТАТОВ
    # =============================================================================
    def process_file(e: ft.FilePickerResultEvent):
        nonlocal result_df
        global last_result_path
        if e.files:
            for f in e.files:
                # Считываем Excel
                df = pd.read_excel(f.path)
                show_preview_excel(df)

                results = []
                for _, row in df.iterrows():
                    # Собираем ионы в моль/л
                    ions = {}
                    if "Cl⁻, мг/л" in row: ions['Cl'] = conc_molar(row["Cl⁻, мг/л"], Mr['Cl'])
                    if "SO₄²⁻ мг/л" in row: ions['SO4'] = conc_molar(row["SO₄²⁻ мг/л"], Mr['SO4'])
                    if "HCO₃, мг/л" in row: ions['HCO3'] = conc_molar(row["HCO₃, мг/л"], Mr['HCO3'])
                    if "Ca²⁺, мг/л" in row: ions['Ca'] = conc_molar(row["Ca²⁺, мг/л"], Mr['Ca'])
                    if "Mg²⁺, мг/л" in row: ions['Mg'] = conc_molar(row["Mg²⁺, мг/л"], Mr['Mg'])
                    if "Na⁺ и K⁺, мг/л" in row: ions['NaK'] = conc_molar(row["Na⁺ и K⁺, мг/л"], Mr['NaK'])
                    if "Ba²⁺, мг/л" in row: ions['Ba'] = conc_molar(row["Ba²⁺, мг/л"], Mr['Ba'])
                    if "Sr²⁺, мг/л" in row: ions['Sr'] = conc_molar(row["Sr²⁺, мг/л"], Mr['Sr'])

                    # Параметры термобарических условий и pH
                    T = row.get("Температура, °С", 20)
                    P = row.get("Давление, атм", 1)
                    pH = row.get("рН", 7.0)

                    # Параметры скважины
                    Q = row.get("Дебит жидкости, м³/сут", 0)
                    water_c = row.get("Обводненность, %", 0)
                    GF = row.get("Газовый фактор, м³/м³", 0)

                    # Параметры нефти и CO2
                    dens_oil = row.get("Плотность нефти, кг/м³", 0)
                    Mroil = row.get("Средняя молекулярная масса пл. нефти, г/моль", 0)
                    y_CO2 = row.get("Мольная доля CO₂ в газовой фазе", 0)
                    X_CO2 = row.get("Мольная доля CO₂ в пл. нефти", 0)

                    # Рассчет масс осадков
                    masses, si_res = calculate_masses(ions, T, P, pH, GF, Q, water_c, y_CO2, X_CO2, Mroil, dens_oil)

                    # Добавляем в результаты с префиксами для разделения
                    result_row = row.to_dict()
                    for salt, mass in masses.items():
                        result_row[f"{salt}, г/л"] = mass
                    for salt, si in si_res.items():
                        # Не добавляем колонку SI для общей строки с суммой масс
                        if salt == 'Общая масса солей':
                            continue
                        result_row[f"SI_{salt}"] = si
                    results.append(result_row)

                # Сохраняем результаты
                result_df = pd.DataFrame(results)
                last_result_path = os.path.join(DOWNLOADS_DIR, f"Результат_{f.name}")
                result_df.to_excel(last_result_path, index=False)
                
                calculate_excel_btn.disabled = False

                # Защита от нажатия, если нет результатов
                if result_df is not None:
                    calculate_excel_btn.on_click = lambda _, df=result_df: show_result_preview(df)
                else:
                    calculate_excel_btn.on_click = None
                
                calculate_excel_btn.style=ft.ButtonStyle(
                    color=ft.Colors.WHITE,
                    bgcolor=ft.Colors.GREEN
                )
                page.update()

    file_picker.on_result = process_file

    # =============================================================================
    # ФУНКЦИЯ ДЛЯ ОТОБРАЖЕНИЯ ПРЕДПРОСМОТРА РЕЗУЛЬТАТОВ В ВИДЕ ТАБЛИЦЫ
    # =============================================================================
    def show_result_preview(df_results):
        nonlocal preview_container

        # Выбираем колонки: №, Скважина, Masses и SI values
        columns_to_show = ['№', 'Скважина']
        
        # Добавляем колонки масс
        mass_columns = [col for col in df_results.columns if col.endswith(', г/л')]
        columns_to_show.extend(mass_columns)
        
        # Добавляем колонки SI
        si_columns = [col for col in df_results.columns if col.startswith('SI_')]
        columns_to_show.extend(si_columns)
        
        # Фильтруем DataFrame
        df_display = df_results[columns_to_show]

        # Формируем DataTable
        columns = [ft.DataColumn(ft.Container(content=ft.Text(col, size=14), alignment=ft.alignment.center, expand=True)) for col in df_display.columns]
        rows = []
        for _, row in df_display.iterrows():
            row_cells = []
            for col in df_display.columns:
                try:
                    val = row[col]
                    if pd.isna(val):
                        display = ''
                    else:
                        # If numeric zero, display as integer zero
                        if isinstance(val, float) and val == 0.0:
                            display = '0'
                        elif col in ['№', 'Скважина']:
                            try:
                                display = str(int(val))
                            except Exception:
                                display = str(val)
                        else:
                            display = str(val)
                except Exception:
                    display = str(row.get(col, ''))
                # Center the cell content
                row_cells.append(ft.DataCell(ft.Container(content=ft.Text(display, size=14), alignment=ft.alignment.center, expand=True)))
            rows.append(ft.DataRow(cells=row_cells))

        table = ft.DataTable(
            columns=columns,
            rows=rows,
            border=ft.border.all(1),
            border_radius=4,
            horizontal_lines=ft.border.BorderSide(1, ft.Colors.BLACK),
            vertical_lines=ft.border.BorderSide(1, ft.Colors.BLACK),
            column_spacing=15,
        )

        preview_container.content = ft.Container(
            ft.Column([
                ft.Text('Результаты расчета', size=22, weight=ft.FontWeight.BOLD),
                ft.Row([table], scroll=ft.ScrollMode.ALWAYS)
            ], horizontal_alignment=ft.CrossAxisAlignment.START, scroll=ft.ScrollMode.ALWAYS)
        )

        # Сбрасываем оформление контейнера
        preview_container.border = None
        preview_container.border_radius = 0
        preview_container.alignment = None
        preview_container.width = None
        preview_container.height = None

        page.update()
        page.open(ft.SnackBar(ft.Text(f"* Ваш файл появился в папке \"Загрузки\""), bgcolor=ft.Colors.GREEN, duration=20000))

    # ============================================================================
    # ФУНКЦИЯ ДЛЯ ОТОБРАЖЕНИЯ ПРЕДПРОСМОТРА ИМПОРТИРОВАННОГО EXCEL
    # ============================================================================
    def show_preview_excel(df: pd.DataFrame):
        # Берем первые 5 строк
        num_rows = max(1, 5)
        num_cols = max(11, len(df.columns))

        df_preview = df.iloc[:num_rows, :num_cols]
        
        # Создаем колонки и строки для DataTable
        columns = [ft.DataColumn(ft.Text(col)) for col in df_preview.columns]
        rows = []
        for _, row in df_preview.iterrows():
            row_cells = []
            for col in df_preview.columns:
                try:
                    if col == 'Дата отбора':
                        val = row[col]
                        if pd.isna(val):
                            cell_text = ''
                        else:
                            try:
                                ts = pd.to_datetime(val)
                                cell_text = ts.strftime('%d.%m.%Y')
                            except Exception:
                                cell_text = str(val)
                    else:
                        val = row[col]
                        if pd.isna(val):
                            cell_text = ''
                        else:
                            # Попытка преобразовать числовые значения в округленные целые для предпросмотра
                            try:
                                cell_text = str(int(round(float(val), 0)))
                            except Exception:
                                cell_text = str(val)
                except Exception:
                    cell_text = str(row.get(col, ''))
                row_cells.append(ft.DataCell(ft.Text(cell_text)))
            rows.append(ft.DataRow(cells=row_cells))
        
        table = ft.DataTable(
            columns=columns,
            rows=rows,
            border=ft.border.all(1),
            border_radius=4,
            horizontal_lines=ft.border.BorderSide(1, ft.Colors.BLACK),
            vertical_lines=ft.border.BorderSide(1, ft.Colors.BLACK),
            column_spacing=20,
        )
        
        table_width = max(1000, 140 * len(df_preview.columns))
        preview_container.width = max(2000, table_width - 100)
        preview_container.content = ft.Container(
            ft.Column([
                ft.Text("Окно предпросмотра данных", size=18, weight=ft.FontWeight.BOLD),
                ft.Row([ft.Container(content=table, padding=5, width=table_width)], scroll=ft.ScrollMode.ALWAYS, wrap=False)
            ], ),
            padding=20
        )
        page.update()

    # =============================================================================
    # СТРАНИЦА EXCEL РЕЖИМА
    # =============================================================================
    def settings_page():
        return ft.Container(
            content=ft.Column(
                [   ft.Row([
                        ft.Text("Режим Excel — для расчётов по нескольким скважинам", size=24, weight=ft.FontWeight.BOLD),
                        theme_btn],
                        alignment=ft.MainAxisAlignment.SPACE_BETWEEN
                        ),
                    ft.Text("Скачайте шаблон Excel, заполните свои данные и загрузите его обратно — программа автоматически всё посчитает и создаст файл с результатами.", size=16),
                    ft.Row(
                        [   
                            ft.ElevatedButton(
                                content=ft.Row([
                                                ft.Text("1. Скачать шаблон", size=16),
                                                ft.Icon(ft.Icons.DOWNLOAD)
                                                ],
                                alignment=ft.MainAxisAlignment.CENTER, spacing=5),
                                on_click=create_template,
                                width=280,
                                height=50,
                                style=ft.ButtonStyle(
                                    color=ft.Colors.WHITE,
                                    bgcolor=ft.Colors.BLUE
                                                    )),
                            
                            ft.ElevatedButton(
                                content=ft.Row([
                                                ft.Text("2. Импорт шаблона", size=16),
                                                ft.Icon(ft.Icons.FOLDER_OPEN)
                                                ],
                                alignment=ft.MainAxisAlignment.CENTER, spacing=5),
                                on_click=lambda _: file_picker.pick_files(allow_multiple=True),
                                width=280,
                                height=50,
                                style=ft.ButtonStyle(
                                    color=ft.Colors.WHITE,
                                    bgcolor=ft.Colors.BLUE
                                                    )
                             ),
                            
                            calculate_excel_btn,

                            ft.ElevatedButton(
                                content=ft.Row([
                                                ft.Text("Очистить", size=16),
                                                ft.Icon(ft.Icons.DELETE)
                                                ],
                                alignment=ft.MainAxisAlignment.CENTER, spacing=5),
                                on_click=lambda e: clear_all(e),
                                width=280,
                                height=50,
                                style=ft.ButtonStyle(
                                    color=ft.Colors.WHITE,
                                    bgcolor=ft.Colors.BLUE
                                                    )
                            )
                        ],
                        spacing=20,
                    ),
                    preview_container
                ],
                spacing=20,
                alignment=ft.MainAxisAlignment.CENTER,
                horizontal_alignment=ft.CrossAxisAlignment.START
            ),
            alignment=ft.alignment.center,
            expand=True,
            padding=10
        )
    
    # =========================================================================
    # СБРОС КОНТЕЙНЕРА ПРЕДПРОСМОТРА
    # =========================================================================
    def clear_all(e):
        preview_container.content = ft.Text(
            "Здесь появится таблица c результатами расчета",
            size=20,
            color=ft.Colors.GREY
        )
        preview_container.border = ft.border.all(1, ft.Colors.GREY)
        preview_container.border_radius = 30
        preview_container.alignment = ft.alignment.center
        preview_container.width = 1200
        preview_container.height = 400

        calculate_excel_btn.disabled = True
        calculate_excel_btn.style=ft.ButtonStyle(
                    color=ft.Colors.WHITE,
                    bgcolor=ft.Colors.GREY
                )

        results_column.controls.clear()

        page.update()
    # =============================================================================
    # --- Основной контейнер для контента страниц ---
    content = ft.Container(expand=True)

    # --- Обработчик переключения вкладок ---
    def on_tab_change(e):
        selected = e.control.selected_index
        if selected == 0:
            content.content = home_page()
        elif selected == 1:
            content.content = settings_page()
        page.update()
    
    # --- Панель навигации ---
    nav_bar = ft.CupertinoNavigationBar(
        bgcolor=ft.Colors.BLUE_GREY_900,
        inactive_color=ft.Colors.GREY_600,
        active_color=ft.Colors.BLUE_400,
        on_change=on_tab_change,
        destinations=[
            ft.NavigationBarDestination(icon=ft.Icons.EDIT_OUTLINED, label="Ручной ввод"),
            ft.NavigationBarDestination(icon=ft.Icons.DESCRIPTION_OUTLINED, label="Excel-режим"),
        ]
    )

    page.navigation_bar = nav_bar

    # --- Первоначальная страница ---
    content.content = home_page()

    # --- Добавляем всё на страницу ---
    page.add(
        ft.Column(
            controls=[content],
            alignment=ft.MainAxisAlignment.SPACE_BETWEEN,
            expand=True
        )
    )

ft.app(target=main)
