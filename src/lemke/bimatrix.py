# bimatrix class

import fractions
import random  # random.seed
import sys

import click
import numpy as np

from . import columnprint, lemke, randomstart, utils


# for debugging
def printglobals():
    globs = [x for x in globals().keys() if not "__" in x]
    for var in globs:
        value = str(globals()[var])
        if not "<" in value:
            print("    " + str(var) + "=", value)


# file format:
# <m> <n>
# m*n entries of A, separated by blanks / newlines
# m*n entries of B, separated by blanks / newlines
#
# blank lines or lines starting with "#" are ignored

# defaults
gamefilename = "game"
gz0 = False
LHstring = ""  # empty means LH not called
seed = -1
trace = -1  # negative: no tracing
accuracy = 1000


# used for both A and B
class payoffmatrix:
    # create matrix from any numerical matrix
    def __init__(self, A):
        AA = np.array(A)
        m, n = AA.shape
        self.numrows = m
        self.numcolumns = n
        self.matrix = np.zeros((m, n), dtype=fractions.Fraction)
        for i in range(m):
            for j in range(n):
                self.matrix[i][j] = utils.tofraction(AA[i][j])
        self.fullmaxmin()

    def __str__(self):
        buf = columnprint.columnprint(self.numcolumns)
        for i in range(self.numrows):
            for j in range(self.numcolumns):
                buf.sprint(str(self.matrix[i][j]))
        out = str(buf)
        out += "\n# max= " + str(self.max) + ", min= " + str(self.min)
        out += ", negshift= " + str(self.negshift)
        return out

    def updatemaxmin(self, fromrow, fromcol):
        m = self.numrows
        n = self.numcolumns
        for i in range(fromrow, m):
            for j in range(fromcol, n):
                elt = self.matrix[i][j]
                self.max = max(self.max, elt)
                self.min = min(self.min, elt)
        self.negshift = int(self.max) + 1
        self.negmatrix = np.full((m, n), self.negshift, dtype=int) - self.matrix

    def fullmaxmin(self):
        self.max = self.matrix[0][0]
        self.min = self.matrix[0][0]
        self.updatemaxmin(0, 0)

    def addrow(self, row):
        self.matrix = np.vstack([self.matrix, row])
        self.numrows += 1
        self.updatemaxmin(self.numrows - 1, 0)

    def addcolumn(self, col):
        self.matrix = np.column_stack([self.matrix, col])
        self.numcolumns += 1
        self.updatemaxmin(0, self.numcolumns - 1)


class bimatrix:
    def __init__(self, filename):
        lines = utils.stripcomments(filename)
        words = utils.towords(lines)
        m = int(words[0])
        n = int(words[1])
        needfracs = 2 * m * n
        if len(words) != needfracs + 2:
            print("in bimatrix file " + repr(filename) + ":")
            print(
                "m=", m, ", n=", n, ", need", needfracs, "payoffs, got", len(words) - 2
            )
            exit(1)
        k = 2
        C = utils.tomatrix(m, n, words, k)
        self.A = payoffmatrix(C)
        k += m * n
        C = utils.tomatrix(m, n, words, k)
        self.B = payoffmatrix(C)

    def __str__(self):
        out = "# m,n= \n" + str(self.A.numrows)
        out += " " + str(self.A.numcolumns)
        out += "\n# A= \n" + str(self.A)
        out += "\n# B= \n" + str(self.B)
        return out

    def createLCP(self):
        m = self.A.numrows
        n = self.A.numcolumns
        lcpdim = m + n + 2
        lcp = lemke.lcp(lcpdim)
        lcp.q[lcpdim - 2] = -1
        lcp.q[lcpdim - 1] = -1
        for i in range(m):
            lcp.M[lcpdim - 2][i] = 1
            lcp.M[i][lcpdim - 2] = -1
        for j in range(m, m + n):
            lcp.M[lcpdim - 1][j] = 1
            lcp.M[j][lcpdim - 1] = -1
        for i in range(m):
            for j in range(n):
                lcp.M[i][j + m] = self.A.negmatrix[i][j]
        for j in range(n):
            for i in range(m):
                lcp.M[j + m][i] = self.B.negmatrix[i][j]
        for i in range(lcpdim):
            lcp.d[i] = 1
        return lcp

    def runLH(self, droppedlabel):
        lcp = self.createLCP()
        lcp.d[droppedlabel - 1] = 0
        tabl = lemke.tableau(lcp)
        tabl.runlemke(silent=True)
        return tuple(getequil(tabl))

    def LH(self, LHstring):
        if LHstring == "":
            return
        m = self.A.numrows
        n = self.A.numcolumns
        lhset = {}
        labels = rangesplit(LHstring, m + n)
        for k in labels:
            eq = self.runLH(k)
            if eq in lhset:
                lhset[eq].append(k)
            else:
                print("label", k, "found eq", str_eq(eq, m, n))
                lhset[eq] = [k]
        print("-------- equilibria found: --------")
        for eq in lhset:
            print(str_eq(eq, m, n), "found by labels", str(lhset[eq]))
        return lhset

    def runtrace(self, xprior, yprior):
        lcp = self.createLCP()
        Ay = self.A.negmatrix @ yprior
        xB = xprior @ self.B.negmatrix
        lcp.d = np.hstack((Ay, xB, [1, 1]))
        tabl = lemke.tableau(lcp)
        tabl.runlemke(silent=True)
        return tuple(getequil(tabl))

    def tracing(self, trace):
        if trace < 0:
            return
        m = self.A.numrows
        n = self.A.numcolumns
        trset = {}
        if trace == 0:
            xprior = uniform(m)
            yprior = uniform(n)
            eq = self.runtrace(xprior, yprior)
            trset[eq] = 1
            trace = 1
        else:
            for k in range(trace):
                if seed >= 0:
                    random.seed(10 * trace * seed + k)
                x = randomstart.randInSimplex(m)
                xprior = randomstart.roundArray(x, accuracy)
                y = randomstart.randInSimplex(n)
                yprior = randomstart.roundArray(y, accuracy)
                eq = self.runtrace(xprior, yprior)
                if eq in trset:
                    trset[eq] += 1
                else:
                    print("found eq", str_eq(eq, m, n), "index", self.eqindex(eq, m, n))
                    trset[eq] = 1
        print("-------- statistics of equilibria found: --------")
        for eq in trset:
            print(trset[eq], "times found ", str_eq(eq, m, n))
        print(trace, "total priors,", len(trset), "equilibria found")

    def eqindex(self, eq, m, n):
        rowset, colset = supports(eq, m, n)
        k, l = len(rowset), len(colset)
        if k != l:
            return 0
        A1 = submatrix(self.A.negmatrix, rowset, colset)
        DA = np.linalg.det(A1)
        B1 = submatrix(self.B.negmatrix, rowset, colset)
        DB = np.linalg.det(B1)
        sign = 2 * (k % 2) - 1
        if DA * DB == 0:
            return 0
        if DA * DB > 0:
            return sign
        return -sign


def uniform(n):
    return np.array([fractions.Fraction(1, n) for j in range(n)])


def getequil(tabl):
    tabl.createsol()
    return tabl.solution[1 : tabl.n - 1]


def str_eq(eq, m, n):
    x = "(" + ",".join([str(x) for x in eq[0:m]]) + ")"
    y = "(" + ",".join([str(x) for x in eq[m : m + n]]) + ")"
    rowset, colset = supports(eq, m, n)
    return x + "," + y + "\n    supports: " + str(rowset) + str(colset)


def supports(eq, m, n):
    rowset = [i for i in range(m) if eq[i] != 0]
    colset = [j for j in range(n) if eq[m + j] != 0]
    return rowset, colset


def submatrix(A, rowset, colset):
    k, l = len(rowset), len(colset)
    B = np.zeros((k, l))
    for i in range(k):
        for j in range(l):
            B[i][j] = A[rowset[i]][colset[j]]
    return B


def rangesplit(s, endrange=50):
    result = []
    for part in s.split(","):
        if part != "":
            if "-" in part:
                a, b = part.split("-")
                a = int(a)
                b = endrange if b == "" else int(b)
            else:
                a = int(part)
                b = a
            a = max(a, 1)
            b = min(b, endrange)
            result.extend(range(a, b + 1))
    return result


# Sentinel: distinguishes "option not passed" from "option passed without value"
_UNSET = object()


@click.command()
@click.argument("gamefilename", default="game")
@click.option(
    "-LH",
    "lh_range",
    default=None,
    metavar="RANGE",
    help="Lemke-Howson with missing labels, e.g. '1,3-5,7-' (omit value = all labels).",
)
@click.option(
    "-trace",
    "trace",
    default=None,
    type=int,
    metavar="NUM",
    help="Tracing procedure; NUM = number of priors (0 = centroid).",
)
@click.option(
    "-seed", "seed_val", default=None, type=int, metavar="NUM", help="Random seed."
)
@click.option(
    "-accuracy",
    "accuracy_val",
    default=1000,
    type=int,
    metavar="N",
    help="Accuracy prior denominator (default: 1000).",
)
@click.option(
    "-decimals",
    "decimals_val",
    default=None,
    type=int,
    metavar="D",
    help="Allowed payoff digits after decimal point (default: 4).",
)
@click.option(
    "-z0", "gz0_flag", is_flag=True, default=False, help="Use z0 entering variable."
)
def main(gamefilename, lh_range, trace, seed_val, accuracy_val, decimals_val, gz0_flag):
    global gz0, LHstring, seed, accuracy

    # -LH with no argument in original code defaulted to "1-" (all labels).
    # With Click, omitting -LH entirely means lh_range is None (LH not called).
    # Passing -LH with an explicit range string works as before.
    # To replicate "-LH with no argument" pass -LH "1-".
    LHstring = lh_range if lh_range is not None else ""

    # -trace with no argument in original code defaulted to 0 (centroid).
    # With Click, omitting -trace means trace stays -1 (tracing not called).
    # Passing -trace 0 explicitly replicates the original "-trace" with no arg.
    trace_val = trace if trace is not None else -1

    # -seed with no argument in original code defaulted to 0.
    # With Click, omitting -seed means seed stays -1 (no fixed seed).
    seed = seed_val if seed_val is not None else -1

    accuracy = accuracy_val
    gz0 = gz0_flag

    if decimals_val is not None:
        utils.setdecimals(decimals_val)

    printglobals()

    G = bimatrix(gamefilename)
    print(G)
    G.LH(LHstring)
    G.tracing(trace_val)


if __name__ == "__main__":
    main()
