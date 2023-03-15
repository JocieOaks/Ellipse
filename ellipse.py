"""
Module for calculating the elliptical projection of a circle onto the viewing plane, used for drawing programs to draw circles in perspective.

Can be run as a script to take a quaternion and a diameter to get the angle and minor axis of a circle with the given diameter being rotated by the given quaternion.

new_quaternion - Create a new unit quaternion with the given axis and angle.
hamilton - Calculate the Hamilton product of two quaternions.
conjugate - Calculate the conjugate of a quaternion.
rotation - Rotate a vector by a unit quaternion.
axis_rotation - Perform two rotations, where the second axis of rotation is also being rotated.
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

def new_quaternion(vector, theta):
    """
    Create a new unit quaternion with the given axis and angle.

    Keyword Arguments:
        vector (float[3]) -- a vector pointing along the axis of rotation.
        theta (float) -- the angle of rotation around the axis in radians.
    Returns:
        float[4] -- a unit quaternion with the given rotation.
    """
    r_2 = vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2
    if(not m.isclose(r_2,1)):
        r = m.sqrt(r_2)
        vector = [vector[i] / r for i in range(0,3)]

    sin = m.sin(theta)
    cos = m.cos(theta)

    return [cos, sin * vector[0], sin * vector[1], sin * vector[2]]

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

def rotation(quaternion, vector):
    """
    Rotate a vector by a unit quaternion.

    Keyword Arguments:
        quaternion (float[4]) -- the unit quaternion about which vector is being rotated.
        vector (float[3]) -- the vector being rotated.
    Returns:
        float[3] -- vector after being rotated by quaternion.
    """
    rotated_vector = hamilton(hamilton(quaternion, [0, vector[0], vector[1], vector[2]]), conjugate(quaternion))
    #hamilton product always returns a quaternion, but the scaler part will be zero, so it is simply a vector.
    return [rotated_vector[i] for i in range(1,4)]

def axis_rotation(quaternion1, axis, theta):
    """
    Caclculate the unit quaternion representing two rotations, where the axis of rotation is also being rotated.

    This is useful for cases such as a gyroscope, where two rings are rotating inside one another, but as the first ring rotates, the axis of the second ring is also rotated.
    Keyword Arguments:
        quaternion (float[4]) -- the unit quaternion representing the initial rotation.
        axis (float[3]) -- the axis of the second rotation.
        theta (float) -- the angle of rotation in radians for the second rotation.
    Returns:
        float[4] -- the quaternion that combines both rotations.
    """
    new_axis = rotation(quaternion1, axis)
    quaternion2 = new_quaternion(new_axis, theta)

    #successive rotation by two quaternions is equal to rotation by the hamilton product of both quaternions.
    return hamilton(quaternion2, quaternion1)


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
    # find one of the two foci of the ellipse then use that to solve the angle and minor axis length
    focus = rotation(quaternion, [0,0,1])
    
    if (m.isclose(focus[1], 0, abs_tol = 0.00001)):
        theta = 0
    elif(m.isclose(focus[0],0, abs_tol = 0.00001)):
        theta = m.pi /2
    else:
        theta = m.atan(focus[1]/ focus[0])

    if (m.isclose(focus[0] ** 2 + focus[1] ** 2, 1)):
       minor = 0
    else:
        minor = m.sqrt(1 - focus[0] ** 2 - focus[1] ** 2)

    return (theta, minor)

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