"""
Module for calculating the elliptical projection of a circle onto the viewing plane, used for drawing programs to draw circles in perspective.

Can be run as a script to take a quaternion and a diameter to get the angle and minor axis of a circle with the given diameter being rotated by the given quaternion.

hamilton - Calculate the Hamilton product of two quaternions.
conjugate - Calculate the conjugate of a quaternion.
rotation - Rotate one quaternion by a unit quaternion.
get_minor_axis - Calculate the length of the minor axis of a described ellipse.
test_ellipse - Determine the goodness of fit for an ellipse.
elliptical_projection - Find the elliptical projection of a rotated circle.
normalize - Normalize a quaternion.
parse_arguments - Parse arguments when run as a script.
main - Main function when run as a script.
"""

import math as m
import argparse
from array import array

#Constants for the positional values of a quarternion array, correlating to their basis.
#Scaler component of a quaternion.
A = 0
#Vector components of a quaternion.
I = 1       
J = 2
K = 3

#The square root of 2 divided by 2, which is the length of the legs of a right triangle with a hypotnuse of legnth 1. 
RAD2_2 = 0.7071067812

def hamilton(quaternion1, quaternion2):
    """
    Calculate the Hamilton product of two quaternions.

    Note that the Hamilton product is not commutative. 
    Keyword Arguments:
        quaternion1 (float[4]) -- the first quaternion.
        quaternion2 (float[4]) -- the second quaternion.
    Returns:
        float[4] -- (quaternion1 * quaternion2)
    """
    a = quaternion1[A] * quaternion2[A] - quaternion1[I] * quaternion2[I] - quaternion1[J] * quaternion2[J] - quaternion1[K] * quaternion2[K]
    i = quaternion1[A] * quaternion2[I] + quaternion1[I] * quaternion2[A] + quaternion1[J] * quaternion2[K] - quaternion1[K] * quaternion2[J]
    j = quaternion1[A] * quaternion2[J] - quaternion1[I] * quaternion2[K] + quaternion1[J] * quaternion2[A] + quaternion1[K] * quaternion2[I]
    k = quaternion1[A] * quaternion2[K] + quaternion1[I] * quaternion2[J] - quaternion1[J] * quaternion2[I] + quaternion1[K] * quaternion2[A]

    return [a, i, j, k]

def conjugate(quaternion):
    """
    Calculate the conjugate of a quaternion.

    Keyword Arguments:
        quaternion (float[4]) -- the quaternion being conjugated.
    Returns:
        float[4] -- quaternion*
    """
    return [quaternion[A], -quaternion[I], -quaternion[J], -quaternion[K]]

def rotation(quaternion1, quaternion2):
    """
    Calculate a quaternion rotated around a unit quaternion.

    Keyword Arguments:
        quaternion1 -- The unit quaternion about which quaternion2 is being rotated.
        quaternion2 -- The quaternion being rotated.
    Returns:
        float[4] -- quaternion2 rotated around quaternion1      aka     quaternion2' = quaternion1 * quaternion2 * quaternion1*
    """
    return hamilton(hamilton(quaternion1,quaternion2),conjugate(quaternion1))

def get_minor_axis(x1, y1, x2, y2, cos, sin):
    """
    Calculate the length of the minor axis of an ellipse with the given parameters.

    Keyword Arguments:
        x1 (float) -- the x position of the first coordinate on the ellipse.
        y1 (float) -- the y position of the first coordinate on the ellipse.
        x2 (float) -- the x position of the second coordinate on the ellipse.
        y2 (float) -- the y position of the second coordinate on the ellipse.
        cos (float) -- the cosine of the angle of the ellipse.
        sin (float) -- the sine of the angle of the ellipse.
    Returns:
        float -- the square of the half length of the minor axis of the ellipse.
    """

    # Based on the equation x^2 / a^2 + y^2 / b^2 = 1 for an ellipse.
    # b is the half length of the major axis and is presumed 1.
    # a is the half length of the minor axis
    # Because the ellipse is at an angle x is replaced with x': x' = x * cos - y * sin
    # and y is replaced with y': y' = x * sin + y * cos

    # There are two possible solutions to this equation, one where the major axis is 1 and one where the minor axis is 1
    # Given that we want the solution where the minor axis is 1, we use the coordinates that are closest to the origin.

    if(m.fabs(x1 * sin + y1 * cos) < m.fabs(x2 * sin + y2 * cos)):
        denom = 1 - m.pow(x1 * sin + y1 * cos, 2)
        if(not m.isclose(denom, 0, abs_tol = 0.00001)):
            return m.pow(x1 * cos - y1 * sin, 2) / denom
        else:
            return 0
    else:
        denom = 1 - m.pow(x2 * sin + y2 * cos, 2)
        if(not m.isclose(denom, 0, abs_tol = 0.00001)):
            return m.pow(x2 * cos - y2 * sin, 2) / denom
        else:
            return 0


def test_ellipse(x, y, sin, cos, minor2):
    """
    Evalute the goodness of fit of an ellipse using a test coordinate.

    Keyword Arguments:
        x (float) -- the x position of the testing coordinate.
        y (float) -- the y position of the testing coordinate.
        sin (float) -- the sine of the ellipse angle.
        cos (float) -- the cosine of the ellipse angle.
        minor2 (float) -- the square of the half legnth of the ellipse' minor axis.
    Returns:
        float -- the absolute value of the ellipse equation. For a coordinate that is on the ellipse, this value should be zero.
    """
    if(not m.isclose(minor2, 0, abs_tol = 1e-10)):
        return m.fabs(m.pow(x * cos - y * sin, 2) / minor2 + m.pow(x * sin + y * cos, 2) - 1)
    else:
        return m.fabs(1 - x * cos - y * sin)

def elliptical_projection(quaternion):
    """
    Find the elliptical projection of a rotated circle.

    Keyword Arguments:
        quaternion float[4] -- the unit quaternion describing the rotation of the circle.
    Returns:
        (float, float) -- the angle of the ellipse in radians, and the half length of the minor axis / major axis. 
                          The major axis is presumed 1, therefore by multiplying the second value by diameter of the circle you will have the full length of minor axis.

    The program defines the ellipse by it's minor axis and angle, because in most drawing programs such as CSP you can define the length of the major and minor axes,
    and then perform a rotation on the ellipse to get it at the correct angle.
    """
    # First definte 3 points that will be on the ellipse by rotating the 0 degree, 45 degree, and 90 degree points of a circle.
    vec1 = rotation(quaternion, [0,1,0,0])
    vec2 = rotation(quaternion, [0,0,1,0])
    vec3 = rotation(quaternion, [0, RAD2_2, RAD2_2, 0])

    # Precalculate several useful values.
    x1 = vec1[I]
    y1 = vec1[J]

    x1_2 = x1 * x1
    y1_2 = y1 * y1

    x2 = vec2[I]
    y2 = vec2[J]

    x2_2 = x2 * x2
    y2_2 = y2 * y2

    x3 = vec3[I]
    y3 = vec3[J]

    # Returns immediately if the circle is orthogonal to the viewing plane and thus appears as a single line.
    if((m.isclose(x1,x2, abs_tol = 1e-5) and m.isclose(y1,y2, abs_tol = 1e-5)) or (m.isclose(x1,-x2, abs_tol = 1e-5) and m.isclose(y1,-y2, abs_tol = 1e-5))):
        if(m.isclose(y1, 0, abs_tol = 1e-5)):
            return (m.pi / 2,0)
        else:
            return (m.atan(x1/y1), 0)

    # Evaluate a system of quadratic equations using the first two coordinates.
    a = x1_2 + y1_2 - x2_2 - y2_2

    b = x1_2 - y1_2 - x2_2 + y2_2 + 2 * x2_2 * y1_2 - 2 * x1_2 * y2_2

    c = 2 * x1 * y1 - 2 * x2 * y2 - 2 * x1 * x2_2 * y1 - 2 * x1 * y1 * y2_2 + 2 * x1_2 * x2 * y2 + 2 * x2 * y1_2 * y2

    a_ = b * b + c * c
    b_ = 2 * a * b
    c_ = a * a - c * c

    determinant = (b_ * b_) - (4 * a_ * c_)
    if(m.isclose(determinant, 0, abs_tol = 0.000001)):
        determinant = 0

    # Negative determinant is irrational, and a_ that is equal to zero creates a divide by zero error.
    if(determinant >= 0 and not m.isclose(a_, 0)):
        u_pos = (-b_ + m.sqrt(determinant)) / (2 * a_)
        u_neg = (-b_ - m.sqrt(determinant)) / (2 * a_)

        if(u_pos > 1 or u_pos < -1):
            theta_pos = m.atan(y1 / x1)
        else:
            theta_pos = m.acos(u_pos) / 2

        if(u_neg > 1 or u_neg < -1):
            theta_neg = m.atan(y1 / x1)
        else:
            theta_neg = m.acos(u_neg) / 2

        # Because the root of a number can be positive or negative there are multiple solutions
        sin_pos = m.sin(theta_pos)
        cos_pos = m.cos(theta_pos)

        minor2_pos = get_minor_axis(x1, y1, x2, y2, cos_pos, sin_pos)
        minor2_pos_ = get_minor_axis(x1, y1, x2, y2, cos_pos, -sin_pos)

        sin_neg = m.sin(theta_neg)
        cos_neg = m.cos(theta_neg)

        minor2_neg = get_minor_axis(x1, y1, x2, y2, cos_neg, sin_neg)
        minor2_neg_ = get_minor_axis(x1, y1, x2, y2, cos_neg, -sin_neg)

        # Use the third coordinate to find the solution that has the closest fit.
        best = [test_ellipse(x3, y3, sin_pos, cos_pos, minor2_pos), theta_pos, minor2_pos]

        test = test_ellipse(x3, y3, -sin_pos, cos_pos, minor2_pos_)

        if(test < best[0]):
            best = [test, -theta_pos, minor2_pos_]

        test = test_ellipse(x3, y3, sin_neg, cos_neg, minor2_neg)

        if(test < best[0]):
            best = [test, theta_neg, minor2_neg]

        test = test_ellipse(x3, y3, -sin_neg, cos_neg, minor2_neg_)

        if(test < best[0]):
            best = [test, -theta_neg, minor2_neg_]

        return (best[1], m.sqrt(best[2]))

    else:
        return (0, 1)


def normalize(quaternion):
    """
    Normalize a quaternion.

    Keyword Arguments:
        quaternion (float[4]) -- the quaternion being normalized.

    Results:
        float[4] -- a unit quaternion that points in the direction of quaternion
    """
    r_2 = quaternion[A] * quaternion[A] + quaternion[I] * quaternion[I] + quaternion[J] * quaternion[J] + quaternion[K] * quaternion[K]
    r = m.sqrt(r_2)
    return [ quaternion[A] / r, quaternion[I] / r, quaternion[J] / r, quaternion[K] / r ]

def parse_arguments():
    """
    Parse the arguments for the script.

    Results:
        namespace -- populated namespace containing the script arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('quaternion', nargs = 4, type = float, help = "4 number quaternion representing the rotation of the circle", default = [])
    parser.add_argument('diameter', type = float, help = "diameter of the circle being rotated")

    return parser.parse_args()

def main():
    """
    Run the script, finding an elliptical projection of a circle with the given rotation and major axis.
    """
    args = parse_arguments()

    (theta, minor) = elliptical_projection(normalize(array('f', args.quaternion)))

    print("angle: ", theta)
    print("minor axis: ", minor * args.diameter * 2)

if __name__ == '__main__':
    main()