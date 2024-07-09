from itertools import product, permutations, repeat
from sys import argv
from copy import copy, deepcopy
import pprint

class Grid:
    def __init__(self, dims: tuple[int], data):
        """Dims is a tuple of n dimensions, data is a n-deep array of n-lists of true,
        false, and none. If the ith value of the point p is true, then there is a stick in
        the ith dimension with its bottom-left corner at p.
        """
        self.dims = dims
        self.data = data

    def __str__(self):
        if len(self.dims) != 2:
            return "Can only print 2d grids."

        m = self.dims[0]
        n = self.dims[1]

        s = ["·", " "] * m + ["·"]
        lines = []
        for _ in range(n):
            lines.append(copy(s))
            lines.append([" "] * (2 * m + 1))
        lines.append(copy(s))

        for i in range(m + 1):
            for j in range(n + 1):
                if self.data[i][j][0] is True:
                    lines[-(2 * j + 1)][2 * i + 1] = "—"
                if self.data[i][j][1] is True:
                    lines[-(2 * j + 2)][2 * i] = "|"
        return "\n".join("".join(line) for line in lines)


    def get_stick(self, point: tuple[int], dim: int) -> bool:
        pointer = self.data
        for coord in point:
            try:
                pointer = pointer[coord]
            except IndexError:
                return None

        return pointer[dim]

    def incr_point(self, point: tuple[int], dim: int) -> tuple[int]:
        new_point = list(point)
        new_point[dim] += 1
        return tuple(new_point)

    def check_square(self, point: tuple[int], i: int, j: int) -> bool:
        if self.get_stick(point, i) and self.get_stick(point, j):
            if self.get_stick(self.incr_point(point, i), j) != self.get_stick(self.incr_point(point, j), i):
                # print(f"Square at {point} in {i}-{j} dimension failed 3 of 4.")
                return False

        return True

    def check_gravity(self, start: tuple[int], offset_dir: int, gravity_dir: int) -> bool:
        pt = start
        # find the first False or None
        while self.get_stick(pt, offset_dir) is True:
            pt = self.incr_point(pt, gravity_dir)

        while self.get_stick(pt, offset_dir) is False:
            pt = self.incr_point(pt, gravity_dir)
            if self.get_stick(pt, offset_dir) is True:
                # print(f"Gravity in {gravity_dir} dimension failed at {pt} in offset direction.")
                return False

        return True

    def check_all(self) -> bool:
        # find all the points, check all the squares and gravity in every direction. too
        # slow but oh well
        for pt in product(*(range(max(1, dim + 1)) for dim in self.dims)):
            for (i, j) in permutations(range(len(self.dims)), 2):
                if not self.check_square(pt, i, j):
                    return False
                if not self.check_gravity(pt, i, j):
                    return False

        return True

# commented out: enumerating grids in 2d
# def improve_inners(m_opts, n_opts, inners):
#     tups = []
#     for m_opt in m_opts:
#         for n_opt in n_opts:
#             tups.append((m_opt, n_opt))
#
#     new_inners = []
#     for inner in inners:
#         for tup in tups:
#             new_inner = copy(inner)
#             new_inner.append(tup)
#             new_inners.append(new_inner)
#
#     return new_inners
#
# def inner_loop(m, n, outers, m_opts):
#     inners = [[]]
#     for j in range(n):
#         n_opts = (False, True)
#         inners = improve_inners(m_opts, n_opts, inners)
#
#     n_opts = (None,)
#     inners = improve_inners(m_opts, n_opts, inners)
#
#     new_outers = []
#     for outer in outers:
#         for inner in inners:
#             new_outer = copy(outer)
#             new_outer.append(inner)
#             new_outers.append(new_outer)
#
#     return new_outers
#
# def enumerate_datas(m, n):
#     outers = [[]]
#     for i in range(m):
#         m_opts = (False, True)
#         outers = inner_loop(m, n, outers, m_opts)
#
#     m_opts = (None,)
#
#     return inner_loop(m, n, outers, m_opts)
#
#
# def enumerate_grids(n, m):
#     for data in enumerate_datas(n, m):
#         yield Grid((n, m), data)
    



def prepare_grid(dims: tuple[int]):
    ret = []
    for dim in dims[::-1]:
        ret = [deepcopy(ret) for _ in range(dim + 1)]

    return ret

def opts_at_point(point, dims):
    opts = [(False, True) for _ in range(len(dims))]
    for i in range(len(dims)):
        if point[i] == dims[i]:
            opts[i] = (None,)

    return [*product(*opts)]

def enumerate_datas(dims: tuple[int]):
    grids = [prepare_grid(dims)]
    
    for point in product(*(range(dim + 1) for dim in dims)):
        opts = opts_at_point(point, dims)
        new_grids = []
        for opt in opts:
            for grid in grids:
                new_grid = deepcopy(grid)
                pointer = new_grid
                for coord in point:
                    pointer = pointer[coord]

                pointer[:] = opt
                
                new_grids.append(new_grid)

        grids = new_grids

    return grids


def enumerate_grids(dims: tuple[int]):
    for data in enumerate_datas(dims):
        yield Grid(dims, data)


if __name__ == "__main__":
    # g = Grid((1, 1, 1),
    #  [[[[False, True, True], [False, True, None]],
    #    [[False, None, True], [True, None, None]]],
    #   [[[None, False, False], [None, True, None]],
    #    [[None, None, True], [None, None, None]]]])
      
    tot = 0
    cnt = 0
    for grid in enumerate_grids(tuple(int(i) for i in argv[1:])):
        tot += 1
        if grid.check_all():
            # pprint.pp(grid.data)
            # print()
    #         print()
    #         print("--------")
    #         print()
            cnt += 1

    print(f"games checked: {tot}")
    print(f"valid games: {cnt}")
