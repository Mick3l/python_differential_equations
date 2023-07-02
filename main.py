import numpy as np
import matplotlib.pyplot as plt
import sympy


if __name__ == '__main__':
    with open('in.txt', "r") as f:
        # считываем значения параметров
        b, tau = float(f.readline()), float(f.readline())
        n, step = int(f.readline()), float(f.readline())
        # считаем начальные условия
        t0, p0 = float(f.readline()), float(f.readline())

    print(f'стартовая точка t = {t0}')
    print(f'количество разбиений - {n} с шагом {step}')

    # аналитическое решение
    p1 = sympy.Function('p')
    t = sympy.symbols('t')
    b1 = sympy.symbols('b')
    tau1 = sympy.symbols('tau')
    diffeq = sympy.Eq((p1(t).diff(t)), b1 * (tau1 - t) * np.e ** (-t / tau1))
    exp = sympy.dsolve(diffeq, ics={p1(0): 0}).rhs

    def p(t_in):
        try:
            return [float(exp.subs('b', b).subs('tau', tau).subs('t', i)) for i in t_in]
        except TypeError:
            return float(exp.subs('b', b).subs('tau', tau).subs('t', t_in))

    # правая часть дифференциального уравнения
    def f(t, p):
        return b * (tau - t) * np.e ** (-t / tau)

    # Численный явный метод Эйлера
    def Euler(h):
        p = p0
        t = t0

        p_array = np.array([p for i in range(n)])
        t_array = p_array.copy()
        t_array[0] = t

        # Вычисляем значения импульса
        for i in range(1, n):
            pi = p + h * f(t, p)
            p = pi
            t += h
            p_array[i] = p
            t_array[i] = t
        return t_array, p_array

    # Построение графика функции p(t), по данным полученным методом Эйлера
    solution1 = Euler(step)
    plt.plot(solution1[0], solution1[1], label="метод Эйлера")
    plt.title("График зависимости p(t)")
    plt.xlabel("t, с")
    plt.ylabel("p(t), кг*м/c")

    # формируем таблицу по скорости и времени и заносим туда данные
    with open('out1.txt', 'w') as out:
        for i in range(len(solution1[0])):
            out.write(f'{solution1[0][i]} \t {solution1[1][i]} \n')

    # Построение графика функции p(t), полученной аналитически
    t = np.linspace(t0, t0 + n * step, 1200)
    plt.plot(t, p(t), "g", label="Аналитическое решение")
    plt.legend(fontsize=15)
    plt.grid()
    plt.show()

    # Оценка порядка точности явного метода Эйлера методом Ричардсона
    m1 = Euler(step * 2)
    n1 = Euler(step * 4)
    p1 = np.log2(abs((solution1[1][-1] - m1[1][-1]) / (m1[1][-1] - n1[1][-1])))
    print('Порядок точности метода Эйлера : ', p1)

    # метод Рунге-Кутты 4 порядка точности

    def Runge_Kutte(h):
        p = p0
        t = t0
        p_array = np.array([p0 for i in range(n)])
        t_array = p_array.copy()
        t_array[0] = t0
        # Вычисляем значения импульса и добавляем их в список
        for i in range(1, n):
            k1 = f(t, p)
            k2 = f(t + h / 2, p + (h * k1) / 2)
            k3 = f(t + h / 2, p + (h * k2) / 2)
            k4 = f(t + h, p + h * k3)
            pi = p + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            p = pi
            t += h
            p_array[i] = pi
            t_array[i] = t
        return t_array, p_array

    # Построение графика функции p(t), по данным полученным методом Рунге-Кутты
    solution2 = Runge_Kutte(step)
    plt.plot(solution2[0], solution2[1], "r", label="Метод Рунге-Кутты")
    plt.title("Графки зависимости p(t)")
    plt.xlabel("t, с")
    plt.ylabel("p(t), м/c")

    # вывод данных в таблицу импульса от времени
    with open('out2.txt', 'w') as out:
        for i in range(len(solution2[0])):
            out.write(f'{solution2[0][i]} \t {solution2[1][i]} \n')

    # Построение графика аналитического решения
    t = np.linspace(t0, t0 + n * step, 1200)
    plt.plot(t, p(t), "g", label="Аналитическое решение")
    plt.legend(fontsize=15, loc=(1.05, 0.5))
    plt.grid()
    plt.show()

    # Оценка точности метода Рунге-Кутты методом Ричардсона
    m2 = Runge_Kutte(step * 2)
    n2 = Runge_Kutte(step * 4)
    p2 = np.log2(abs((solution2[1][-1] - m2[1][-1]) / (m2[1][-1] - n2[1][-1])))
    print('Порядок точности метода Рунге-Кутты : ', p2)

    # График зависимости численной ошибки от шага для явного метода Эйлера:
    steps = np.array(
        [0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 3, 0, 0.5, 1, 2, 3, 5, 10, 20, 50, 200, 500, 1000, 2000])
    # Массив со значениями ошибки
    error_euler = np.array([Euler(i)[1][-1] for i in steps])
    error_euler -= p(t0 + n * step)
    error_runge = np.array([Runge_Kutte(i)[1][-1] for i in steps])
    error_runge -= p(t0 + n * step)

    # График зависимости ошибки от шага
    plt.plot(steps, error_euler)
    plt.plot(steps, error_runge, "r")
    plt.title("Графки зависимости численной ошибки (Err), от шага (h)")
    plt.xlabel("h")
    plt.ylabel("Err")
    plt.grid()
    plt.show()

    # Метода логарифмических невязок
    log_error_euler = np.log(error_euler)
    log_error_runge = np.log(error_runge)
    log_steps = np.log(steps)
    plt.plot(log_steps, log_error_euler)
    plt.plot(log_steps, log_error_runge, "r")
    plt.title("Определение порядка точности методом логарифмических невязок")
    plt.xlabel("Логарифм шага")
    plt.ylabel("Логарифмическая невязка ")
    plt.grid()
    plt.show()

    # основываясь на графики метода логарифмических невязом наилучший шаг = 0.05
    #  исследуем зависимость решения от параметра tau задачи
    n = 100  # чтобы значения функции были не слищком большие
    taus = [-5, -1, 0.5, 1, 5, 15]
    for i in taus:
        tau = i
        solution2 = Runge_Kutte(step)
        plt.plot(solution2[0], solution2[1])
        plt.title("Графки зависимости p(t)")
        plt.xlabel("t, с")
        plt.ylabel("p(t), м/c")
    plt.grid()
    plt.show()
