"""
Script for visualizing the ellipse module. Draws 3 rotating concentric rings.

print_ellipse - Draw an ellipse in a graphics window.
print_axis - Draw the axis of rotation of a rotating ring.
parse_arguments - Parse the arguments for the script.
main - Run the script.
"""

import ellipse
import graphics as g
import math as m
import time as t
import argparse

def print_ellipse(win, offset, major, minor, theta):
    """
    Draw an ellipse in a graphics window.

    Keyword Arguments:
        win (GraphWin) -- the window in which to draw the ellipse.
        offset (int) -- the number of pixels to shift the origin of the ellipse down and to the right.
        major (float) -- the length of the major axis of the ellipse.
        minor (float) -- the length of the minor axis of the ellipse.
        theta (float) -- the angle of the ellipse' major axis, in radians.
    """

    sin = m.sin(theta)
    cos = m.cos(theta)

    old = g.Point(minor * cos + offset, minor * sin + offset)

    count = 48

    for i in range(1, count):
        x = minor * m.cos(i * m.pi * 2 / count)
        y = major * m.sin(i * m.pi * 2 / count)

        new = g.Point(x * cos - y * sin + offset, x * sin + y * cos + offset)

        line = g.Line(old,new)
        line.draw(win)
        old = new

    new = g.Point(minor * cos + offset, minor * sin + offset)

    line = g.Line(old,new)
    line.draw(win)

def print_axis(vector, win, offset, diameter_inner, dimater_outer):
    """
    Draw the axis of rotation of a rotating ring.

    Keyoword Arguments:
        vector (float[3]) -- the vector along which the axis lies.
        win (GraphWin) -- the window in which to draw the ellipse.
        offset (int) -- the number of pixels to shift the origin of the ellipse down and to the right.
        diameter_inner (float) -- the diameter of the rotating ring.
        diameter_out (float) -- the dimater of the ring within which the ring is rotating.
    """
    p1 = g.Point(vector[0] * diameter_inner + offset, vector[1] * diameter_inner + offset)
    p2 = g.Point(vector[0] * dimater_outer + offset, vector[1] * dimater_outer + offset)

    line = g.Line(p1,p2)
    line.draw(win)

    p1 = g.Point(-vector[0] * diameter_inner + offset, - vector[1] * diameter_inner + offset)
    p2 = g.Point(-vector[0] * dimater_outer + offset, -vector[1] * dimater_outer + offset)

    line = g.Line(p1,p2)
    line.draw(win)

def parse_arguments():
    """
    Parse the arguments for the script.

    Results:
        namespace -- populated namespace containing the script arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('angle1', type = float, help = "angle in degrees of the first ring's rotational axis")
    parser.add_argument('speed1', type = float, help = "speed of the first ring's rotation, defined as rotations pers second")
    parser.add_argument('angle2', type = float, help = "angle in degrees of the second ring's rotational axis")
    parser.add_argument('speed2', type = float, help = "speed of the second ring's rotation, defined as rotations pers second")
    parser.add_argument('angle3', type = float, help = "angle in degrees of the third ring's rotational axis")
    parser.add_argument('speed3', type = float, help = "speed of the third ring's rotation, defined as rotations pers second")

    return parser.parse_args()

def main():
    """
    Display 3 rings rotating around one another until the window is clicked on or closed by the user.
    """
    size = 800

    win = g.GraphWin("Test", size, size)

    args = parse_arguments()

    major1 = 375
    major2 = 300
    major3 = 225

    angular_velocity1 = m.pi * args.speed1 / 20
    angular_velocity2 = m.pi * args.speed2 / 20
    angular_velocity3 = m.pi * args.speed3 / 20

    angle1 = -args.angle1 * m.pi / 180
    angle2 = -args.angle2 * m.pi / 180
    angle3 = -args.angle3 * m.pi / 180

    # These vectors are technically quaternions, but the scaler portion is 0.
    vector1 = [m.cos(angle1), m.sin(angle1), 0]
    vector2 = [m.cos(angle2), m.sin(angle2), 0]
    vector3 = [m.cos(angle3), m.sin(angle3), 0]

    i = 0

    while(win.checkMouse() == None and not win.isClosed()):

        quat1 = ellipse.new_quaternion(vector1, angular_velocity1 * i)
        quat2 = ellipse.axis_rotation(quat1, vector2, angular_velocity2 * i)
        quat3 = ellipse.axis_rotation(quat2, vector3, angular_velocity3 * i)

        measurements1 = ellipse.elliptical_projection(quat1)
        measurements2 = ellipse.elliptical_projection(quat2)
        measurements3 = ellipse.elliptical_projection(quat3)

        win.delete("all")
        print_ellipse(win, size / 2, major1, measurements1[1] * major1, measurements1[0])
        print_ellipse(win, size / 2, major2, measurements2[1] * major2, measurements2[0])
        print_ellipse(win, size / 2, major3, measurements3[1] * major3, measurements3[0])

        print_axis(vector1, win, size/2, major1, major1 + 75)
        print_axis(ellipse.rotation(quat1,vector2), win, size/2, major2, major1)
        print_axis(ellipse.rotation(quat2,vector3), win, size/2, major3, major2)

        i += 1

        t.sleep(0.05)

    win.close()

if __name__ == '__main__':
    main()