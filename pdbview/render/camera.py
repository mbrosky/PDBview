"""3D camera and projection"""

import numpy as np
from typing import Tuple


class Camera:
    """3D camera with rotation, zoom, and pan"""

    def __init__(self, distance: float = 50.0):
        self.distance = distance
        self.rotation_x = 0.0
        self.rotation_y = 0.0
        self.rotation_z = 0.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.zoom = 1.0

    def rotate(self, dx: float, dy: float, dz: float = 0.0):
        """Rotate camera"""
        self.rotation_x += dx
        self.rotation_y += dy
        self.rotation_z += dz

    def pan(self, dx: float, dy: float):
        """Pan camera"""
        self.pan_x += dx
        self.pan_y += dy

    def zoom_in(self, factor: float = 1.1):
        """Zoom in"""
        self.zoom *= factor

    def zoom_out(self, factor: float = 1.1):
        """Zoom out"""
        self.zoom /= factor

    def reset(self):
        """Reset to default view"""
        self.rotation_x = 0.0
        self.rotation_y = 0.0
        self.rotation_z = 0.0
        self.pan_x = 0.0
        self.pan_y = 0.0
        self.zoom = 1.0

    def get_rotation_matrix(self) -> np.ndarray:
        """Get combined rotation matrix"""
        # Rotation around X axis
        rx = np.radians(self.rotation_x)
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(rx), -np.sin(rx)],
            [0, np.sin(rx), np.cos(rx)]
        ])

        # Rotation around Y axis
        ry = np.radians(self.rotation_y)
        Ry = np.array([
            [np.cos(ry), 0, np.sin(ry)],
            [0, 1, 0],
            [-np.sin(ry), 0, np.cos(ry)]
        ])

        # Rotation around Z axis
        rz = np.radians(self.rotation_z)
        Rz = np.array([
            [np.cos(rz), -np.sin(rz), 0],
            [np.sin(rz), np.cos(rz), 0],
            [0, 0, 1]
        ])

        # Combine rotations: Rz * Ry * Rx
        return Rz @ Ry @ Rx

    def project(self, points: np.ndarray, width: int, height: int) -> np.ndarray:
        """
        Project 3D points to 2D screen coordinates.

        Args:
            points: Nx3 array of 3D coordinates
            width: Screen width
            height: Screen height

        Returns:
            Nx2 array of 2D screen coordinates
        """
        if len(points) == 0:
            return np.array([])

        # Apply rotation
        R = self.get_rotation_matrix()
        rotated = points @ R.T

        # Apply zoom and pan
        scale = self.zoom * min(width, height) / (2 * self.distance)

        # Orthographic projection (simple but works)
        x_2d = rotated[:, 0] * scale + width / 2 + self.pan_x
        y_2d = -rotated[:, 1] * scale + height / 2 + self.pan_y  # Flip Y for screen coords
        z_2d = rotated[:, 2]  # Keep Z for depth sorting

        return np.column_stack([x_2d, y_2d, z_2d])

    def project_point(self, point: np.ndarray, width: int, height: int) -> Tuple[float, float, float]:
        """Project a single 3D point to 2D"""
        result = self.project(point.reshape(1, 3), width, height)
        return tuple(result[0])
