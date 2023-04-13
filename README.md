# Ellipse
Small module for calculating the properties of a rotating ring and performing quaternion mathematics in regards to them. The primary motivation of this was for use with art programs, such as Photoshop and Clip Studio Paint, where an ellipse can be defined by it's major and minor axes, and then rotated by some angle in degrees. This allows an artist to make accurate animations of a rotating ring, such as with a gyroscope.

<p align="center">
  <img src="https://github.com/JocieOaks/Ellipse/blob/master/Gyroscope.gif" alt="Example" height="300" width="auto">
</p>

### Prerequisites
The visual.py script requires [Zelle Graphics](https://mcsp.wartburg.edu/zelle/python/graphics.py).

### Usage
The ellipse.py module can be used as a library for performing Quaternion calculations and visualizing a rotating ring, or at can be run as a script taking the quaternion angle that represents the rotation of a ring that was originally parallel to the screen and the diameter, and returns the minor axis and rotation of the resulting ellipse.

### License
This project is under an MIT license.
