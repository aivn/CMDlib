#!/usr/bin/python3
"""
Модуль для работы с dat-данными. 

Под dat-данными принимается следующая структура:
----------
#:x y z    <- заголовочный файл   
1 2 3      |
4 5 6      | <- dat-блок
7 8 9      |
            <- dat-блоки разделяются пустыми строками
7 8 9
10 11 1
3 4 5
----------

'shape' данных -- это размер данных без учета колонок и с учетом размерности.
В примере выше количество колонок равно трём, а 'shape' равен (2, 3).
Заметьте, что размерность 'shape' может равняться только 1 или 2, 
но это не помешает эмулировать массивы больших размерностей (которые
пока не поддерживаются).
"""
import itertools
import numpy as np
import pandas as pd


def parse_line(line):
    """
    Разделяет линию с помощью 'split', переводит все во 'float'.
    Если в строке нет элементов, возвращает None.
    """
    splitted = line.split()
    if not splitted:
        return None

    result = []
    for elem in splitted:
        try:
            elem = float(elem)
        except ValueError:
            pass

        result.append(elem)

    return [[elem] for elem in result]


def try_parse_names(line, start="#:"):
    """
    Пытается запарсить заголовочную линию. Заголовочная линия начинается с 'start'.
    Если линия заголовочная -> возвращает [name0, name1, ...], иначе -> None.
    """
    if line[0:len(start)] == start:
        return [name for name in line[len(start):].split()]
    return None


def parse_block(lines):
    """
    Парсит один dat-блок из итерируемого объекта 'lines'.
    """
    result = []

    try:
        values = parse_line(next(lines))
    except StopIteration:
        return result

    result = values

    while True:
        try:
            values = parse_line(next(lines))
        except StopIteration:
            values = None

        if values is None:
            break

        for i, elem in enumerate(values):
            result[i] += elem

    return result


def parse(lines, start="#:"):
    """
    Парсит итерируемый объект 'lines'. 
    При наличии заголовка -> {name: index}, [[[e00, e01, ...], [e10, e11, ...], ...], ...],
    иначе -> None, [[[e00, e01, ...], [e10, e11, ...], ...], ...].

    #:x y z           
    1 2 3             data[names['x']] = [1, 4, 7]
    4 5 6       ->    data[names['y']] = [2, 5, 6]
    7 8 9             data[names['z']] = [3, 6, 9]

    В случае наличия несколько dat-блоков второй индекс отвечает за номер блока.
    """
    result = []

    try:
        probably_names_line = next(lines)
    except StopIteration:
        return None, np.asarray(result)

    names = try_parse_names(probably_names_line, start)
    if names is None:
        lines = itertools.chain([probably_names_line], lines)
        lines = iter(lines)

    block = parse_block(lines)
    if not block:
        return names, np.asarray(result)
    result = [[elem] for elem in block]

    while True:
        block = parse_block(lines)
        if not block:
            break

        for i, elem in enumerate(block):
            result[i].append(elem)

    number_or_object = True
    for col in result:
        for block in col:
            for elem in block:
                number_or_object = number_or_object and (
                    isinstance(elem, int) or isinstance(elem, float))

    for i, elem in enumerate(result):
        if len(elem) == 1:
            result[i] = elem[0]

    if number_or_object:
        return names, np.asarray(result)
    else:
        return names, np.asarray(result, dtype=object)


def parse_from_input(input, start="#:"):
    """
    см. parse()
    """
    return parse(iter(input.readlines()), start)


def np2dat(filename, names, data, start="#:"):
    """
    Превращает 'names' и 'data' в dat-данные. 
    """
    file = open(filename, "w")
    result = ""
    
    if names is not None:
        file.write(start + " ".join(names) + "\n")
    
    shape = data.shape[1:]
    cols = data.shape[0]
    
    # 1D
    if len(shape) == 1:
        for i in range(shape[-1]):
            for col in range(cols):
                file.write(str(data[col][i]) + " ")
            file.write("\n")
    # 2D
    if len(shape) == 2:
        for j in range(shape[-2]):
            for i in range(shape[1]):
                for col in range(cols):
                   file.write(str(data[col][j, i]) + " ")
                file.write("\n")
            file.write("\n")
    
    if len(shape) > 2:
        raise ValueError(f"не допускается dim больше 2: dim={len(shape)}")
    
    file.close()
