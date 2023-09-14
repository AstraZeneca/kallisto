# test/test_grid.py
from kallisto.grid import getLebedevLaikovGrid
from kallisto.grid import gridSize


def test_user_can_create_grid_sizes():
    for n in range(23):
        size = gridSize[n]
        grid, weights = getLebedevLaikovGrid(n)
        assert len(weights) == size
        assert len(grid[:, 1]) == size
        assert len(grid[1, :]) == 3
