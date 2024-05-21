"""
    Parameters of the algorithms
"""

class CurveFitMethodParameters:
    def __init__(self):
        self.min_curve_length = 10
        self.max_time_gap = 60

        self.iteration = _CurveFitIterationParameters()
        self.grid_search = _CurveFitGridSearchParameters()


class _CurveFitStepParameter:
    def __init__(self, init_step, step_limit, lower_bound=None, upper_bound=None):
        self.init_step = init_step
        self.step_limit = step_limit
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound


class _CurveFitIterationParameters:
    def __init__(self):
        self.lat = _CurveFitStepParameter(100e3, 1)
        self.lon = _CurveFitStepParameter(100e3, 1)
        self.alt = _CurveFitStepParameter(10,    1, 0, 3000)
        self.off = _CurveFitStepParameter(3000,  1)
        self.dft = _CurveFitStepParameter(0.1,   0.001)

        self.iteration_limit = 500
        self.repeats = 3


class _CurveFitGridSearchParameter:
    def __init__(self, steps, grid_size, lower_bound=None, upper_bound=None):
        self.steps = steps
        self.grid_size = grid_size
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound


class _CurveFitGridSearchParameters:
    def __init__(self):
        self.lat = _CurveFitGridSearchParameter([10], 10)
        self.lon = _CurveFitGridSearchParameter([10], 10)
        self.alt = _CurveFitGridSearchParameter([10], 30, 0, 3000)
        self.off = _CurveFitGridSearchParameter([0], 6)
        self.dft = _CurveFitGridSearchParameter([0], 6)

        self.grid_search_depth = len(self.lat.steps)

        # all steps must have the same length
        assert len(self.lat.steps) == len(self.lon.steps) == len(self.alt.steps) == len(self.off.steps) == len(self.dft.steps)
