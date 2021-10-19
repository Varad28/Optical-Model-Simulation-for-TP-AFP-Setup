# Varad Kulkarni
# Reflection Functions:

import numpy as np
import sympy as sp
from sympy import Eq
from sympy.abc import x, y
import matplotlib.pyplot as plt


def coeff_extract(equations):
    coefficients = list()
    for i in equations:
        p = i.as_poly()
        xy = p.coeff_monomial(x), p.coeff_monomial(y)
        coefficients.append(xy)
    return coefficients


def subs_ref_eqns(previous_solutions, previous_equations):
    normal_vector_components = coeff_extract(previous_equations)
    normal_vector_components = np.array(normal_vector_components, dtype=np.float64)
    unit_normal_vector_components = [i / np.linalg.norm(i) for i in normal_vector_components]
    vector_components = [[i[1], -i[0]] for i in unit_normal_vector_components]  # Forming ray vector of previous rays
    normal_vector = np.array([0, 1])  # Normal vector to substrate
    # Calculating reflected ray vector: r = i - 2(i.n)n
    l1 = [np.dot(i, normal_vector) for i in vector_components]
    l2 = [i * 2 for i in l1]
    l3 = [i * normal_vector for i in l2]
    reflected_ray_vectors = [i - j for i, j in zip(vector_components, l3)]
    reflected_ray_normal_vectors = [[-i[1], i[0]] for i in reflected_ray_vectors]  # Normal vectors of new reflected rays
    reflected_ray_equations = sp.Function('reflected_ray_equations')
    x, y = sp.symbols('x y')
    reflected_ray_equations = [Eq(i[0] * (x - j[0]) + i[1] * (y - j[1]), 0) for i, j in zip(reflected_ray_normal_vectors, previous_solutions)]
    return reflected_ray_equations


def roller_ref_eqns(previous_solutions, previous_equations, k):
    normal_vector_components = coeff_extract(previous_equations)
    normal_vector_components = np.array(normal_vector_components, dtype=np.float64)
    unit_normal_vector_components = [i / np.linalg.norm(i) for i in normal_vector_components]
    unit_incident_vectors = np.array([(i[1], -i[0]) for i in unit_normal_vector_components])  # Forming ray vector using normal vector components
    normal_vectors = [[i[0] - 0, i[1] - k] for i in previous_solutions]
    normal_vectors = np.array(normal_vectors, dtype=np.float64)
    unit_normal_vectors = [i / np.linalg.norm(i) for i in normal_vectors]

    l4 = [np.dot(i, j) for i, j in zip(unit_incident_vectors, unit_normal_vectors)]
    l5 = [i * 2 for i in l4]
    l6 = list()
    for i, j in zip(l5, unit_normal_vectors):
        l7 = list()
        for q in j:
            l7.append(i * q)
        l6.append(l7)
    l8 = [i - j for i, j in zip(unit_incident_vectors, l6)]
    reflected_ray_vectors = [list(i) for i in l8]
    reflected_ray_normal_vectors = [[-i[1], i[0]] for i in reflected_ray_vectors]
    reflected_ray_equations = sp.Function('reflected_ray_equations')
    x, y = sp.symbols('x y')
    reflected_ray_equations = [Eq(i[0] * (x - j[0]) + i[1] * (y - j[1]), 0) for i, j in zip(reflected_ray_normal_vectors, previous_solutions)]
    return reflected_ray_equations


def it_ref_eqns(previous_solutions, previous_equations, a2, b2, pq):
    normal_vector_components = coeff_extract(previous_equations)
    normal_vector_components = np.array(normal_vector_components, dtype=np.float64)
    unit_normal_vector_components = [i / np.linalg.norm(i) for i in normal_vector_components]
    unit_incident_vectors = [(i[1], -i[0]) for i in unit_normal_vector_components]  # Forming ray vector using normal vector components
    normal_vectors = np.array([b2 - pq[1], pq[0] - a2])
    normal_vectors = np.array(normal_vectors, dtype=np.float64)
    unit_normal_vectors = normal_vectors / np.linalg.norm(normal_vectors)

    l9 = [np.dot(i, unit_normal_vectors) for i in unit_incident_vectors]
    l10 = [i * 2 for i in l9]
    l11 = [i * unit_normal_vectors for i in l10]
    reflected_ray_vectors = [i - j for i, j in zip(unit_incident_vectors, l11)]
    reflected_ray_normal_vectors = [[-i[1], i[0]] for i in reflected_ray_vectors]
    reflected_ray_equations = sp.Function('reflected_ray_equations')
    x, y = sp.symbols('x y')
    reflected_ray_equations = [Eq(i[0] * (x - j[0]) + i[1] * (y - j[1]), 0) for i, j in zip(reflected_ray_normal_vectors, previous_solutions)]
    return reflected_ray_equations


def to_subs(equations_reflected_from_tape_or_roller):
    substrate = sp.Function('substrate')
    substrate = Eq(y, 0)
    coordinates = list()
    for i in equations_reflected_from_tape_or_roller:
        l12 = sp.solve([i, substrate], [x, y])
        l13 = list(dict.values(l12))
        coordinates.append(l13)
    return coordinates


def subs_to_rol_it(equations_reflected_from_substrate, a2, b2, pq, k):
    tape = sp.Function('tape')
    roller = sp.Function('roller')
    x, y = sp.symbols('x, y')
    tape = Eq((b2 - pq[1]) * (x - a2) + (pq[0] - a2) * (y - b2), 0)
    roller = Eq(x ** 2 + (y - k) ** 2 - k ** 2, 0)
    substrate_to_tape_coordinates = list()
    substrate_to_tape_equations = list()
    substrate_to_roller_equations = list()
    all_roller_solns = dict()
    for i in equations_reflected_from_substrate:
        l18 = sp.solve([i, tape], [x, y])  # Solving rays and tape
        l19 = list(dict.values(l18))
        if l19[0] >= pq[0] and l19[1] >= pq[1]:
            substrate_to_tape_coordinates.append(l19)
            substrate_to_tape_equations.append(i)
        elif l19[0] < pq[0] and l19[1] < pq[1]:
            l20 = sp.solve([i, roller], [x, y])  # Solving rays and roller
            l21 = {i: l20}
            all_roller_solns.update(l21)
    roller_real_solns = list()
    for i in all_roller_solns:
        for j in all_roller_solns[i]:
            if type(j[0]) is not sp.core.add.Add or type(j[1]) is not sp.core.add.Add:
                roller_real_solns.append(j)
        substrate_to_roller_equations.append(i)
    # Filtering required points of reflections on Roller:
    l22 = list()
    for i in roller_real_solns:
        l22.append(i[0])  # Appending all x values in l22
    l23 = iter(l22)
    l24 = list()
    for i in l23:
        rx = max(i, next(l23))  # Comparing every 2 x coordinates and -
        l24.append(rx)  # -storing max in s1m
    l25 = list()
    for i in roller_real_solns:  # Same for y coordinate
        l25.append(i[1])
    l26 = iter(l25)
    l27 = list()
    for j in l26:
        ry = min(j, next(l26))
        l27.append(ry)
    substrate_to_roller_coordinates = [[i, j] for [i, j] in zip(l24, l27)]
    return substrate_to_tape_equations, substrate_to_tape_coordinates, substrate_to_roller_equations, substrate_to_roller_coordinates


# Plotting Functions:


def plt_incidence(solution_coordinates, equations, a0, b0):
    coefficients = coeff_extract(equations)
    for i, j in zip(solution_coordinates, coefficients):  # Plotting rays incident on incoming tape
        i = list(np.float_(i))
        x1, y1 = np.meshgrid(np.linspace(i[0], a0, 50), np.linspace(i[1], b0, 50))
        plt.contour(x1, y1, j[0] * (x1 - a0) + j[1] * (y1 - b0), [0], linewidths=0.5, colors='r')


def plt_to_subs(solution_coordinates, previous_solutions, equations):
    coefficients = coeff_extract(equations)
    end_points = [one + two for one, two in zip(previous_solutions, solution_coordinates)]
    for i, j in zip(end_points, coefficients):
        i = list(np.float_(i))
        x1, y1 = np.meshgrid(np.linspace(i[0], i[2], 50), np.linspace(i[3], i[1], 50))
        plt.contour(x1, y1, j[0] * (x1 - i[0]) + j[1] * (y1 - i[1]), [0], linewidths=0.5, colors='lightsalmon')


def plt_from_subs(solution_coordinates, equations):
    coefficients = coeff_extract(equations)
    previous_solutions = to_subs(equations)
    end_points = [one + two for one, two in zip(previous_solutions, solution_coordinates)]
    for i, j in zip(end_points, coefficients):
        if i[0] > i[2]:
            a = i[2]
            b = i[0]
        elif i[2] > i[0]:
            a = i[0]
            b = i[2]
        if i[1] > i[3]:
            c = i[3]
            d = i[1]
        elif i[3] > i[1]:
            c = i[1]
            d = i[3]
        a = (np.float_(a))
        b = (np.float_(b))
        c = (np.float_(c))
        d = (np.float_(d))
        x1, y1 = np.meshgrid(np.linspace(a, b, 50), np.linspace(c, d, 50))
        plt.contour(x1, y1, j[0] * (x1 - i[2]) + j[1] * (y1 - i[3]), [0], linewidths=0.5, colors='lightsalmon')
        a = 0.0
        b = 0.0
        c = 0.0
        d = 0.0
