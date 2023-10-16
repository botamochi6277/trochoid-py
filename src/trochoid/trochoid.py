#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file is for draw outer/inner trochoids
"""
import math
import numpy as np


def path_length(x, y):
    """
    Return length of a path.

    Parameters
    ----------
    x : array_like
        x-position of points of the path.
    y : array_like
        y-position of points of the path.
        Array size of x and y shouled be equaled.

    Returns
    -------
    s : float
        length of the path.

    """
    dx = np.diff(x)
    dy = np.diff(y)
    ds = dx * dx + dy * dy
    ds = np.sqrt(ds)
    s = np.sum(ds)
    return s


def lcm(x, y):
    """
    Return least common multiple.

    https://note.nkmk.me/python-gcd-lcm/

    Parameters
    ----------
    x : int
        integer.
    y : int
        annother integer.

    Returns
    -------
    z : int
        least common multiple of x and y.

    """
    return (x * y) // math.gcd(x, y)


def polygon(n, r, cx=0.0, cy=0.0, orient=0, ax=None, *args, **kwargs):
    """
    Calculate vertex positions of a regular polygon.

    Parameters
    ----------
    n : int
        the number of vertexs of the polygon
    r : double
        radius of a circumcircle of the polygon
    cx :
        center x-position of the circumcircle
    cy :
        center y-position of the circumcircle
    orient: double
        orientation of a polygon

    Returns
    -------
    x : array_like
        vertex x-position
    y : array_like
        vertex y-position
    """
    theta = np.linspace(0 + orient, 2 * np.pi + orient, num=n + 1)
    x = r * np.cos(theta) + cx
    y = r * np.sin(theta) + cy
    return (x, y)


def ctrochoid(rc, rm, rd, rmax=None, outer=True, num=50,
              orient=0, *args, **kwargs):
    d = 100
    i = 0
    x0 = 0
    y0 = 0
    x = np.empty(num)
    y = np.empty(num)

    while d > rc / 100:
        theta = np.linspace(2 * np.pi * i + orient, 2 *
                            np.pi * (i + 1) + orient, num=num)
        tmp_x = np.empty(num)
        tmp_y = np.empty(num)
        if outer is True:
            tmp_x = (rc + rm) * np.cos(theta) - rd * \
                np.cos(((rc + rm) / rm) * theta)
            tmp_y = (rc + rm) * np.sin(theta) - rd * \
                np.sin(((rc + rm) / rm) * theta)
        else:
            tmp_x = (rc - rm) * np.cos(theta) + rd * \
                np.cos(((rc - rm) / rm) * theta)
            tmp_y = (rc - rm) * np.sin(theta) - rd * \
                np.sin(((rc - rm) / rm) * theta)

        if i == 0:
            x = tmp_x
            y = tmp_y
            x0 = x[0]
            y0 = y[0]
        else:
            x = np.append(x, tmp_x)
            y = np.append(y, tmp_y)

        d = np.linalg.norm(np.array([x[-1] - x0, y[-1] - y0]), 2)
        i = i + 1
        if i > 100:
            print('ERROR!')
            break
    if rmax is not None:
        if outer is True:
            m = rc + rm + rd
            x = rmax / m * x
            y = rmax / m * y
        else:
            if rc > rm:
                m = rc - rm + rd
                x = rmax / m * x
                y = rmax / m * y
            else:
                m = -rc + rm + rd
                x = rmax / m * x
                y = rmax / m * y

    return (x, y)


def trochoid(px, py, rm, rd, right=True, rmax=None,
             orient=0, *args, **kwargs):
    """
    Calculate trochoid on an optional path.

    Parameters
    ----------
    px : array_like
        x-position of points on a path.
    py : array_like
        y-position of points on a path.
        Array size of x and y shouled be equaled.
    rm : double
        radius of a rolling circle.
    rd : duble
        radius of a drawing circle.
    right: boolean
        when right is true, the rolling circle is on the right of the path.
        when right is false, the rolling circle is on the left of the path.
    orient: double
        initial angle of the rolling circle

    Returns
    -------
    x : array_like
        x-position of points on a drawing path.
    y : array_like
        y-position of points on a drawing path.
    """
    x = np.zeros(len(px))
    y = np.zeros(len(py))
    s = 0  # total rolling length.
    theta = orient  # angle of the rolling circle
    n0 = np.array([[1], [0]])

    rot = np.pi / 2
    if right is True:
        rot = -rot
    r_mat = np.matrix(
        [[np.cos(rot), -np.sin(rot)],
         [np.sin(rot), np.cos(rot)], ]
    )

    for i in range(len(px)):
        ds = 0  # delta-s, partial rolling length.
        if i > 0:
            ds = np.linalg.norm(
                np.array([px[i] - px[i - 1], py[i] - py[i - 1]]))

        d_theta = ds / rm  # partial rolling angle.

        s = s + ds

        theta = theta + d_theta

        # t : tangental vector on the path
        if (i - 1) < 0:
            t = np.array([px[i + 1] - px[i], py[i + 1] - py[i]])
        elif (i + 1) >= len(px):
            t = np.array([px[i] - px[i - 1], py[i] - py[i - 1]])
        else:
            t = np.array([px[i + 1] - px[i - 1], py[i + 1] - py[i - 1]]) * 0.5

        t = t / (np.linalg.norm(t) + 1e-9)  # normalize
        n = np.dot(r_mat, np.reshape(t, (2, 1)))  # normal vector on the path
        if i == 0:
            n0 = -n
        # position of the center of the rolling circle.
        pm = np.array([[px[i]], [py[i]]]) + rm * n

        r_ort = np.matrix(
            [[np.cos(theta), -np.sin(theta)],
             [np.sin(theta), np.cos(theta)], ]
        )

        # position of the drawing point
        p_d = pm + rd * np.dot(r_ort, n0)
        x[i] = p_d.item(0)
        y[i] = p_d.item(1)

    return (x, y)


# if __name__ == '__main__':
#     import matplotlib.pyplot as plt
#     plt.style.use('seaborn-colorblind')
#     plt.style.use('seaborn-whitegrid')
#     plt.rcParams['figure.figsize'] = 600 / 72, 600 / 72
#     plt.rcParams["font.size"] = 16

#     fig = plt.figure()

#     # Hypocycloid & Epicycloid
#     fig.add_subplot(2, 2, 1)
#     x, y = polygon(64, 1)
#     fig.axes[0].plot(x, y)
#     x, y = ctrochoid(1, 1 / 6.0, 1 / 6.0, num=400)
#     fig.axes[0].plot(x, y)
#     x, y = ctrochoid(1, 1 / 6.0, 1 / 6.0, outer=False, num=400)
#     fig.axes[0].plot(x, y)
#     fig.axes[0].set_title('Hypocycloid & Epicycloid')
#     fig.axes[0].set_aspect('equal', 'box')

#     # Hypotrochoid & Epitrochoid
#     fig.add_subplot(2, 2, 2)
#     x, y = polygon(64, 1)
#     fig.axes[1].plot(x, y)
#     x, y = ctrochoid(1, 1 / 6.0, 1.5 / 6.0, num=400)
#     fig.axes[1].plot(x, y)
#     x, y = ctrochoid(1, 1 / 6.0, 1.5 / 6.0, outer=False, num=400)
#     fig.axes[1].plot(x, y)
#     fig.axes[1].set_title('Hypotrochoid & Epitrochoid')
#     fig.axes[1].set_aspect('equal', 'box')

#     # Trochoid on Epitrochoid
#     fig.add_subplot(2, 2, 3)
#     x0, y0 = polygon(64, 1)
#     fig.axes[2].plot(x0, y0)
#     x1, y1 = ctrochoid(rc=1, rm=1 / 6, rd=0.8 * 1 /
#                        6, num=1024, outer=False)
#     fig.axes[2].plot(x1, y1)
#     m = 6
#     n = 12 / m  # num. rotation of the rolling circle.
#     rm = path_length(x1, y1) / (2 * np.pi) / n
#     x, y = trochoid(px=np.tile(x1[::-1], m),
#                     py=np.tile(y1[::-1], m), rm=rm, rd=1.8 * rm)
#     fig.axes[2].plot(x, y)

#     fig.axes[2].set_title('Trochoid on Epitrochoid')
#     fig.axes[2].set_aspect('equal', 'box')

#     # Trochoid on Sine Curve
#     fig.add_subplot(2, 2, 4)
#     x = np.linspace(-1, 1, num=100)
#     y = np.sin(x)
#     fig.axes[3].plot(x, y)
#     xt, yt = trochoid(px=x, py=y, rm=0.1, rd=0.2)
#     fig.axes[3].plot(xt, yt)
#     fig.axes[3].set_title('Trochoid on Sine Curve')

#     plt.show()
