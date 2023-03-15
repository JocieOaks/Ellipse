import ellipse
import graphics as g
import math as m
import time as t
import argparse

def printellipse(win, offset, major, minor, theta):
    """
    Print an ellipse in a graphics window.

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

def anchors(quaternion, win, offset, diameter):
    """
    Print the three anchor points from which elliptical_projection calculates the ellipse.

    Keyword Arguments:
        quaternion (float[4]) -- the quaternion that the circle has been rotated by to create the ellipse.
        win (GraphWin) -- the window in which to draw the anchor points.
        offset (int) -- the number of pixels to shift the origin of the ellipse down and to the right.
        diameter (float) -- the diameter of the circle, which is also the major axis of the ellipse.
    """
    vec1 = ellipse.rotation(quaternion, [0,1,0,0])
    vec2 = ellipse.rotation(quaternion, [0,0,1,0])
    vec3 = ellipse.rotation(quaternion, [0, ellipse.RAD2_2, ellipse.RAD2_2, 0])

    x1 = vec1[1]
    y1 = vec1[2]

    g.Circle(g.Point(x1 * diameter + offset, -y1 * diameter + offset), vec1[3] * 5 + 20).draw(win)

    x2 = vec2[1]
    y2 = vec2[2]
    g.Circle(g.Point(x2 * diameter + offset, -y2 * diameter + offset), vec2[3] * 5 + 20).draw(win)

    x3 = vec3[1]
    y3 = vec3[2]

    g.Circle(g.Point(x3 * diameter + offset, -y3 * diameter + offset), vec3[3] * 5 + 20).draw(win)

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

    speed1 = m.pi * args.speed1 / 20
    speed2 = m.pi * args.speed2 / 20
    speed3 = m.pi * args.speed3 / 20

    angle1 = args.angle1 * m.pi / 180
    angle2 = args.angle2 * m.pi / 180
    angle3 = args.angle3 * m.pi / 180

    # These vectors are technically quaternions, but the scaler portion is 0.
    vector1 = [0, m.cos(angle1), m.sin(angle1), 0]
    vector2 = [0, m.cos(angle2), m.sin(angle2), 0]
    vector3 = [0, m.cos(angle3), m.sin(angle3), 0]

    i = 0

    while(win.checkMouse() == None and not win.isClosed()):
        sin = m.sin(speed1 * i)
        cos = m.cos(speed1 * i)
        quat1 = [cos, sin * vector1[1], sin * vector1[2], sin * vector1[3]]

        vector = ellipse.rotation(quat1,vector2)

        sin = m.sin(speed2 * i)
        cos = m.cos(speed2 * i)

        quat2 = [cos, sin * vector[1], sin * vector[2] , sin * vector[3]]

        vector = ellipse.rotation(ellipse.hamilton(quat2,quat1),vector3)

        sin = m.sin(speed3 * i)
        cos = m.cos(speed3 * i)

        quat3 = [cos, sin * vector[1], sin * vector[2] , sin * vector[3]]

        win.delete("all")
        measurements1 = ellipse.elliptical_projection(quat1)

        anchors(quat1, win, size/2, major1)

        quat2_1 = ellipse.hamilton(quat2,quat1)
        measurements2 = ellipse.elliptical_projection(quat2_1)
        anchors(quat2_1, win, size/2, major2)

        quat3_2_1 = ellipse.hamilton(quat3, quat2_1)
        measurements3 = ellipse.elliptical_projection(quat3_2_1)
        anchors(quat3_2_1,win, size/2, major3)



        printellipse(win, size / 2, major1, measurements1[1] * major1, measurements1[0])
        printellipse(win, size / 2, major2, measurements2[1] * major2, measurements2[0])
        printellipse(win, size / 2, major3, measurements3[1] * major3, measurements3[0])

        i += 1

        t.sleep(0.05)

    win.close()

if __name__ == '__main__':
    main()