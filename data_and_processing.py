import numpy as np
import pandas as pd
import math

import time

# Класс хранения данных
# Методов обработки и получения данных
class DataAndProcessing:
    def __init__(self):
        # Константы
        self.D = 0.0103
        self.a = 0.382
        self.a_6_power = self.a ** 6  # Для расчета сил
        self.constant_before_sum = 12 * self.D * self.a_6_power  # Для расчета сил
        self.particle_size = self.a
        # для скоростей
        self.k = 1.380649  # * 10**(-23) # Постоянная Больцмана [Дж/К]
        self.mass = 0.03995  # кг/моль

        self.step_x = self.a  # Шаг расстановки частиц по x
        self.step_y = self.a * math.sqrt(3) / 2  # Шаг расстановки частиц по y

        # Границы расчетной ячейки
        self.indent_by_x = self.step_x / 2
        self.indent_by_y = self.step_y / 2

        self.start_x = None
        self.start_of_border_x = None
        self.end_of_border_x = None

        self.start_y = None
        self.start_of_border_y = None
        self.end_of_border_y = None

        # Количеств частиц и их координаты
        self.number_elements_on_side = None
        self.temperature = None

        # Частицы, их координаты, скорости и т.д. Индексы - номер частицы.
        title_particles = [
            "x",  # Координата по x
            "y",  # Координата по y

            "Vx",  # Скорость по x
            "Vy",  # Скорость по y

            "Fx",  # Сила по x
            "Fy",  # Сила по y
        ]
        self.particles = pd.DataFrame(columns=title_particles)

        self.energy_kinetic_mass = []

    # (*) Чистит все данные (альтернатива:  self.data = self.data.head(0))
    def clear_data(self):
        self.particles = self.particles.head(0)

    # ОБЩИЙ АЛГОРИТМ
    # Инициализирует начальное состояние
    def initial_state(
            self,
            number_elements_on_side: int,
            start_x=0, start_y=0,
            temperature=300  # Кельвинов
    ):
        self.number_elements_on_side = number_elements_on_side
        self.start_x = start_x
        self.start_y = start_y

        # Расставляем частицы. Находим координаты каждой
        self.coordinate_all_x(start_x, self.step_x, number_elements_on_side)
        self.coordinate_all_y(start_y, self.step_y, number_elements_on_side)

        # Находим границы расчетной ячейки
        self.calculation_boundaries_calculation_cell()

        # Находим скорости
        self.temperature = temperature
        self.speed_by_coordinates()

        # Считаем силы
        self.force(self.particles)

        # Считаем кинетическую энергию
        self.energy_kinetic()
        print(self.energy_kinetic_mass)

    # Итерационный процесс расчета
    def state_calculation(self, time_step=4, number_of_steps=20):
        time_step = 10 ** (-time_step)
        start_time = time.perf_counter()
        for i in range(number_of_steps):
            self.particles = self.step_calculation(self.particles, time_step)
            self.energy_kinetic()
            if i%5==0:
                self.speed_normalization()
        print(self.energy_kinetic_mass)
        end_time = time.perf_counter()
        print(end_time - start_time, "seconds")


    # МЕТОДЫ ПРОЦЕССА ОБРАБОТКИ
    # (0) Расчет границ расчетной ячейки
    def calculation_boundaries_calculation_cell(self):
        # Границы расчетной ячейки
        self.start_of_border_x = self.start_x - self.indent_by_x
        self.end_of_border_x = self.start_x + (self.number_elements_on_side - 1) * self.step_x + self.indent_by_x

        self.start_of_border_y = self.start_y - self.indent_by_y
        self.end_of_border_y = self.start_y + (self.number_elements_on_side - 1) * self.step_y + self.indent_by_y

    # (1) Расстановка частиц - получение координат
    # по x, возвращая массив координат всех частиц
    def coordinate_all_x(self, coordinate_start, step, number_elements_on_side):
        # Создаем пустой массив(столбец) выходных данных
        all_coordinate = np.empty(0, dtype="float64")
        # Подготавливаем четные и не четные строки
        line_even = np.arange(number_elements_on_side, dtype="float64") * step + coordinate_start
        line_odd = np.arange(number_elements_on_side - 1, dtype="float64") * step + step / 2 + coordinate_start
        # По строчно проходим, к итоговому результату прибавляем строки,
        # для новой четной-четная линия и для не четной - не четная
        for line in range(number_elements_on_side):
            if line % 2 == 0:
                all_coordinate = np.append(all_coordinate, line_even)
            else:
                all_coordinate = np.append(all_coordinate, line_odd)
        self.particles.x = all_coordinate

    # по y, возвращая массив координат всех частиц
    def coordinate_all_y(self, coordinate_start, step, number_elements_on_side):
        # Создаем пустой массив(столбец) выходных данных
        all_coordinate = np.empty(0, dtype="float64")

        # По строчно проходим
        for line in range(number_elements_on_side):
            # Количество элементов в строке. Четное - все, не четное - минус один
            amount_elements = 0
            if line % 2 == 0:
                amount_elements = number_elements_on_side
            else:
                amount_elements = number_elements_on_side - 1

            # Добавляем новую строку к итоговому результату
            all_coordinate = np.append(all_coordinate,
                                       # Создаем сторку, выбранного количества элементов
                                       # С одинаковой координатой
                                       np.zeros(amount_elements, dtype="float64") + (line * step + coordinate_start))

        self.particles.y = all_coordinate

    # (2) Задание начальных скоростей
    def speed_by_coordinates(self):
        particle_number = self.particles.index.size
        # Задаем модуль начальной скорости

        """
        mV^2   d
        ---- = - * k * Temperature; d - степень свободы
          2    2

        V = sqrt( d*k*T / m )  
        """
        v_start = math.sqrt((2 * self.k * self.temperature) / self.mass)

        # Модуль начальной скорости x, как случайно равномерно распределенные числа[0,1) * 2pi
        self.particles.Vx = np.random.uniform(0, 1, particle_number) * 2 * np.pi
        # Модуль y = x
        self.particles.Vy = self.particles.Vx
        # Берем cos для x; sin для y и перемножаем с амплитудой - модуль начальной скорости
        self.particles.Vx = v_start * np.cos(self.particles.Vx)
        self.particles.Vy = v_start * np.sin(self.particles.Vy)

        # Нормировка по скорости
        self.speed_normalization()

    # (2*) Нормировка скоростей
    def speed_normalization(self):
        # Нормировка по импульсу
        particle_number = self.particles.index.size
        # Находим средние скорости проекций
        average_speed_x = np.sum(self.particles.Vx) / particle_number
        average_speed_y = np.sum(self.particles.Vy) / particle_number
        # Вычитаем их всех проекций среднее
        self.particles.Vx = self.particles.Vx - average_speed_x
        self.particles.Vy = self.particles.Vy - average_speed_y

    # (3) Расчет сил
    # (3.1) Расчет для одной строки
    def line_force(self, i):
        # (1) Разница xi со всеми x
        difference_x = self.particles.x[i] - self.particles.x
        difference_y = self.particles.y[i] - self.particles.y

        # (2) Квадрат радиуса i - элемента к всем другим
        radius_2_power = difference_x * difference_x + difference_y * difference_y

        # (3) Сила i, для каждой координаты
        # start_time = time.perf_counter ()
        radius_6_power = radius_2_power ** 3
        const_weight = (self.a_6_power / (radius_6_power) - 1) / (radius_2_power * radius_6_power)
        force_x = np.nansum((const_weight) * (difference_x))
        force_y = np.nansum((const_weight) * (difference_y))
        return pd.DataFrame(data={
            "Fx": force_x,
            "Fy": force_y,
        }, index=[i])

    # (3.1) Расчет всех сил
    def force(self, particles):
        particle_number = particles.index.size
        force_x = np.zeros(particle_number, dtype="float64")
        force_y = np.zeros(particle_number, dtype="float64")
        # Перебираем i индекс, по xi yi
        for i in range(particle_number):
            # (1) Разница xi со всеми x
            difference_x = particles.x[i] - particles.x
            difference_y = particles.y[i] - particles.y

            # (2) Квадрат радиуса i - элемента к всем другим
            radius_2_power = difference_x * difference_x + difference_y * difference_y

            # (3) Сила i, для каждой координаты
            radius_6_power = radius_2_power ** 3
            const_weight = (self.a_6_power / radius_6_power - 1) / (radius_2_power * radius_6_power)
            force_x[i] = np.nansum(const_weight * difference_x)
            force_y[i] = np.nansum(const_weight * difference_y)

        particles["Fx"] = self.constant_before_sum * force_x
        particles["Fy"] = self.constant_before_sum * force_y

    # Расчет следующих координат
    # Один шаг итерации
    def step_calculation(self, particles, time_step):

        new_particles = pd.DataFrame(columns=particles.columns)

        new_particles.x = particles.x + \
                          particles.Vx * time_step + \
                          particles.Fx * ((time_step * time_step) / (2 * self.mass))

        new_particles.y = particles.y + \
                          particles.Vy * time_step + \
                          particles.Fy * ((time_step * time_step) / (2 * self.mass))

        DataAndProcessing.check_bounds(self.start_of_border_x, self.end_of_border_x,
                                       new_particles.x)
        DataAndProcessing.check_bounds(self.start_of_border_y, self.end_of_border_y,
                                       new_particles.y)

        self.force(new_particles)

        new_particles.Vx = particles.Vx + (new_particles.Fx + particles.Fx) * (time_step / (2 * self.mass))
        new_particles.Vy = particles.Vy + (new_particles.Fy + particles.Fy) * (time_step / (2 * self.mass))


        return new_particles

    # Проверка выхода частицы за границы
    @staticmethod
    def check_bounds(x_start, x_end, data):
        period = x_end - x_start
        print(period)
        data[data < x_start] = data[data < x_start] + period
        data[data > x_end] = data[data > x_end] - period

    # Считает кинетическую энергию
    def energy_kinetic(self):
        self.energy_kinetic_mass.append(
            np.nansum(self.particles.Vx * self.particles.Vx + self.particles.Vy * self.particles.Vy) * self.mass / 2)