# Varad Kulkarni
# Thermal Model Input Functions:

import itertools
import numpy as np
import sympy as sp
import pandas as pd


# Iterating for range of discretised elements:
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


# Calculating Heat Flux on Substrate:
def calc_subs_heat_flux(element_length, final_value, solution_coordinates, number_of_rays):
    substrate_length = round(final_value)
    n_elements = round(substrate_length / element_length)
    ranges_subs = np.linspace(0, final_value, n_elements)  # All values of nodes on discretised element (Substrate)
    des_sub = [[i, j] for i, j in pairwise(ranges_subs)]  # Coordinates of nodes on discretised element (Substrate)
    heat_flux_list = list()
    x1 = list()
    x2 = list()
    for j in des_sub:
        counter = 0
        heat_flux = 0
        for i in solution_coordinates:
            if j[0] <= i[0] < j[1]:
                counter += 1
                heat_flux = (2.4 / (element_length ** 2 * number_of_rays)) * counter
        heat_flux_list.append(heat_flux)
        x1.append(j[0])
        x2.append(j[1])
    return heat_flux_list, x1, x2,


# Calculating Heat Flux on Incoming Tape:
def calc_it_heat_flux(element_length, solution_coordinates, a2, b2, pq, number_of_rays):
    tape_length = round(sp.sqrt((a2 - pq[0]) ** 2 + (b2 - pq[1]) ** 2))
    n_elements = round(tape_length / element_length)  # No.of elements the IT has to be divided into for the given element length
    ranges_itx = np.linspace(pq[0], a2, n_elements)  # All x values of nodes on discretised element (Substrate)
    ranges_ity = np.linspace(pq[1], b2, n_elements)  # All y values of nodes on discretised element (Substrate)
    des_itx = [[i, j] for i, j in pairwise(ranges_itx)]
    des_ity = [[i, j] for i, j in pairwise(ranges_ity)]
    des_it = [index + des_ity[i] for i, index in enumerate(des_itx)]  # (x1,x2,y1,y2)
    heat_flux_list = list()
    x1 = list()
    y1 = list()
    x2 = list()
    y2 = list()
    for i in des_it:
        counter = 0
        heat_flux = 0
        for j in solution_coordinates:
            if i[0] <= j[0] < i[1] and i[2] <= j[1] < i[3]:
                counter += 1
                heat_flux = (2.4 / (element_length ** 2 * number_of_rays)) * counter
        heat_flux_list.append(heat_flux)
        x1.append(i[0])
        y1.append(i[2])
        x2.append(i[1])
        y2.append(i[3])
    return heat_flux_list, x1, y1, x2, y2


# Calculating Heat Flux on Roller:
def calc_roller_heat_flux(element_length, solution_coordinates, k, arc_length, arc_angle, number_of_rays):
    n_elements = round(arc_length / element_length)
    n_elements = int(n_elements)
    arc_angle = float(arc_angle)
    angle = np.linspace(0, arc_angle, n_elements)
    ranges_rolx = [0 + (k * sp.sin(i)) for i in angle]
    ranges_roly = [k - (k * sp.cos(i)) for i in angle]
    des_rolx = [[i, j] for i, j in pairwise(ranges_rolx)]
    des_roly = [[i, j] for i, j in pairwise(ranges_roly)]
    des_rol = [i + des_roly[index] for index, i in enumerate(des_rolx)]  # (x1,x2,y1,y2)
    heat_flux_list = list()
    x1 = list()
    y1 = list()
    x2 = list()
    y2 = list()
    for i in des_rol:
        counter = 0
        heat_flux = 0
        for j in solution_coordinates:
            if i[0] <= j[0] < i[1] and i[2] <= j[1] < i[3]:
                counter += 1
                heat_flux = (2.4 / (element_length ** 2 * number_of_rays)) * counter
        heat_flux_list.append(heat_flux)
        x1.append(i[0])
        y1.append(i[2])
        x2.append(i[1])
        y2.append(i[3])
    return heat_flux_list, x1, y1, x2, y2


# Exporting Heat Flux data to Excel:
def export_heat_flux(data_frame_name, x1, y1, x2, y2, heat_flux, file_name: str):
    data_frame_name = pd.DataFrame({
        'X1': x1,
        'Y1': y1,
        'X2': x2,
        'Y2': y2,
        'Heat Flux': heat_flux
    })
    data_frame_name.to_csv(file_name, index=False, header=False)
