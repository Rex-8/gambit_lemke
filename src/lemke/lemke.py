# LCP solver

import fractions
import math  # gcd
import sys

import click

from . import columnprint, utils

# global defaults
lcpfilename = "lcp"
outfile = lcpfilename + ".out"
filehandle = sys.stdout
verbose = False
silent = False
z0 = False


# LCP data M,q,d
class lcp:
    # create LCP either with given n or from file
    def __init__(self, arg):
        if isinstance(arg, int):  # arg is an integer
            n = self.n = arg
            self.M = [[]] * n
            for i in range(n):
                self.M[i] = [0] * n
            self.q = [0] * n
            self.d = [0] * n
        else:  # assume arg is a string = name of lcp file
            # create LCP from file
            filename = arg
            lines = utils.stripcomments(filename)
            # flatten into words
            words = utils.towords(lines)
            if words[0] != "n=":
                printout(
                    "lcp file",
                    repr(filename),
                    "must start with 'n=' lcpdim, e.g. 'n= 5', not",
                    repr(words[0]),
                )
                exit(1)
            n = int(words[1])
            self.n = n
            self.M = [[]] * n
            for i in range(n):
                self.M[i] = [0] * n
            self.q = [0] * n
            self.d = [0] * n
            needfracs = n * n + 2 * n
            if len(words) != needfracs + 5:
                printout("in lcp file " + repr(filename) + ":")
                printout(
                    "n=",
                    n,
                    ", need keywords 'M=' 'q=' 'd=' and n*n + n + n =",
                    needfracs,
                    "fractions, got",
                    len(words) - 5,
                )
                exit(1)
            k = 2  # index in words
            while k < len(words):
                if words[k] == "M=":
                    k += 1
                    self.M = utils.tomatrix(n, n, words, k)
                    k += n * n
                elif words[k] == "q=":
                    k += 1
                    self.q = utils.tovector(n, words, k)
                    k += n
                elif words[k] == "d=":
                    k += 1
                    self.d = utils.tovector(n, words, k)
                    k += n
                else:
                    printout("in lcp file " + repr(filename) + ":")
                    printout("expected one of 'M=' 'q=' 'd=', got", repr(words[k]))
                    exit(1)
            return

    def __str__(self):
        n = self.n
        M = self.M
        q = self.q
        d = self.d
        m = columnprint.columnprint(n)
        m.makeLeft(0)
        m.sprint("M=")
        m.newline()
        for i in range(n):
            for j in range(n):
                m.sprint(str(M[i][j]))
        m.sprint("q=")
        m.newline()
        for i in range(n):
            m.sprint(str(q[i]))
        m.sprint("d=")
        m.newline()
        for i in range(n):
            m.sprint(str(d[i]))
        return "n= " + str(n) + "\n" + str(m)

    #######  end of class lcp


class tableau:
    # filling the tableau from the LCP instance Mqd
    def __init__(self, Mqd):
        self.n = Mqd.n
        n = self.n
        self.scalefactor = [0] * (n + 2)  # 0 for z0, n+1 for RHS
        self.A = [[]] * n
        for i in range(n):
            self.A[i] = [0] * (n + 2)
        self.determinant = 1
        self.lextested = [0] * (n + 1)
        self.lexcomparisons = [0] * (n + 1)
        self.pivotcount = 0
        self.solution = [fractions.Fraction(0)] * (2 * n + 1)  # all vars
        self.bascobas = [0] * (2 * n + 1)
        self.whichvar = [0] * (2 * n + 1)
        for i in range(n + 1):  # variables Z(i) all cobasic
            self.bascobas[i] = n + i
            self.whichvar[n + i] = i
        for i in range(n):  # variables W(i+1) all basic
            self.bascobas[n + 1 + i] = i
            self.whichvar[i] = n + 1 + i
        # determine scale factors, lcm of denominators
        for j in range(n + 2):
            factor = 1
            for i in range(n):
                if j == 0:
                    den = Mqd.d[i].denominator
                elif j == n + 1:  # RHS
                    den = Mqd.q[i].denominator
                else:
                    den = Mqd.M[i][j - 1].denominator
                # least common multiple
                factor *= den // math.gcd(factor, den)
            self.scalefactor[j] = factor
            # fill in column j of A
            for i in range(n):
                if j == 0:
                    den = Mqd.d[i].denominator
                    num = Mqd.d[i].numerator
                elif j == n + 1:  # RHS
                    den = Mqd.q[i].denominator
                    num = Mqd.q[i].numerator
                else:
                    den = Mqd.M[i][j - 1].denominator
                    num = Mqd.M[i][j - 1].numerator
                self.A[i][j] = (factor // den) * num
            self.determinant = -1
        return

    def __str__(self):
        out = "Determinant: " + str(self.determinant)
        n = self.n
        tabl = columnprint.columnprint(n + 3)
        tabl.makeLeft(0)
        tabl.sprint("var")  # headers
        for j in range(n + 1):
            tabl.sprint(self.vartoa(self.whichvar[j + n]))
        tabl.sprint("RHS")
        tabl.sprint("scfa")  # scale factors
        for j in range(n + 2):
            if j == n + 1:  # RHS
                tabl.sprint(str(self.scalefactor[n + 1]))
            elif self.whichvar[j + n] > n:  # col  j  is some  W
                tabl.sprint("1")
            else:
                tabl.sprint(str(self.scalefactor[self.whichvar[j + n]]))
        tabl.newline()  # blank line
        for i in range(n):
            tabl.sprint(self.vartoa(self.whichvar[i]))
            for j in range(n + 2):
                s = str(self.A[i][j])
                if s == "0":
                    s = "."  # replace 0 by dot
                tabl.sprint(s)
        out += "\n" + str(tabl)
        out += "\n" + "-----------------end of tableau-----------------"
        return out

    def vartoa(self, v):  # variable as as string w1..wn or z0..zn
        if v > self.n:
            return "w" + str(v - self.n)
        else:
            return "z" + str(v)

    def createsol(self):  # get solution from current tableau
        n = self.n
        for i in range(2 * n + 1):
            row = self.bascobas[i]
            if row < n:  # i is a basic variable
                num = self.A[row][n + 1]
                if i <= n:  # computing Z(i)
                    num *= self.scalefactor[i]
                self.solution[i] = fractions.Fraction(
                    num, self.determinant * self.scalefactor[n + 1]
                )
            else:  # i is nonbasic
                self.solution[i] = fractions.Fraction(0)

    def outsol(self):  # string giving solution, after createsol()
        n = self.n
        sol = columnprint.columnprint(n + 2)
        sol.sprint("basis=")
        for i in range(n + 1):
            if self.bascobas[i] < n:  #  Z(i) is a basic variable
                s = self.vartoa(i)
            elif i > 0 and self.bascobas[n + i] < n:  #  W(i) is a basic variable
                s = self.vartoa(n + i)
            else:
                s = "  "
            sol.sprint(s)
        sol.sprint("z=")
        for i in range(2 * n + 1):
            sol.sprint(str(self.solution[i]))
            if i == n:  # new line since printouting slack vars  w  next
                sol.sprint("w=")
                sol.sprint("")  # no W(0)
        return str(sol)

    def assertbasic(self, v, info):  # assert that v is basic
        if self.bascobas[v] >= self.n:
            printout(info, "Cobasic variable", self.vartoa(v), "should be basic")
            exit(1)
        return

    def assertcobasic(self, v, info):  # assert that v is cobasic
        if self.bascobas[v] < self.n:
            printout(info, "Cobasic variable", self.vartoa(v), "should be cobasic")
            exit(1)
        return

    def docupivot(self, leave, enter):  # leave, enter in VARS
        self.assertbasic(leave, "docupivot")
        self.assertcobasic(enter, "docupivot")
        s = "leaving: " + self.vartoa(leave).ljust(5)
        s += "entering: " + self.vartoa(enter)
        printout(s)
        return

    def raytermination(self, enter):
        printout("Ray termination when trying to enter", self.vartoa(enter))
        printout(self)
        printout("Current basis not an LCP solution:")
        self.createsol()
        printout(self.outsol())
        exit(1)

    def testtablvars(self):  # msg only if error, continue
        n = self.n
        for i in range(2 * n + 1):
            if self.bascobas[self.whichvar[i]] != i:
                for j in range(2 * n + 1):
                    if j == i:
                        printout("First problem for j=", j, ":")
                    printout(
                        f"j={j} self.bascobas[j]={self.bascobas[j]} self.whichvar[j]={self.whichvar[j]}"
                    )
                break
        return

    def complement(self, v):  # Z(i),W(i) are complements
        n = self.n
        if v == 0:
            printout("Attempt to find complement of z0")
            exit(1)
        if v > n:
            return v - n
        else:
            return v + n

    def outstatistics(self):
        n = self.n
        lext = self.lextested
        stats = columnprint.columnprint(n + 2)
        stats.makeLeft(0)
        stats.sprint("lex-column")
        for i in range(n + 1):
            stats.iprint(i)
        stats.sprint("times tested")
        for i in range(n + 1):
            stats.iprint(lext[i])
        if lext[0] > 0:  # otherwise never a degeneracy
            stats.sprint("% of pivots")
            for i in range(0, n + 1):
                stats.iprint(round(lext[i] * 100 / self.pivotcount))
            stats.sprint("avg comparisons")
            for i in range(n + 1):
                if lext[i] > 0:
                    x = round(self.lexcomparisons[i] * 10 / lext[0])
                    stats.sprint(str(x / 10.0))
                else:
                    stats.sprint("-")
        printout(stats)

    def lexminvar(self, enter):
        n = self.n
        A = self.A
        self.assertcobasic(enter, "Lexminvar")
        col = self.bascobas[enter] - n  # entering tableau column
        leavecand = []  # candidates(=rows) for leaving var
        for i in range(n):  # start with positives in entering col
            if A[i][col] > 0:
                leavecand.append(i)
        if leavecand == []:
            self.raytermination(enter)
        if len(leavecand) == 1:
            z0leave = self.bascobas[0] == leavecand[0]
        j = 0  # going through j = 0..n
        while len(leavecand) > 1:
            if j > n:
                printout("lex-minratio test failed")
                exit(1)
            self.lextested[j] += 1
            self.lexcomparisons[j] += len(leavecand)
            if j == 0:
                testcol = n + 1  # RHS
            else:
                testcol = self.bascobas[n + j] - n  # tabl col of W(j)
            if testcol != col:
                if testcol >= 0:
                    newcand = [leavecand[0]]
                    for i in range(1, len(leavecand)):
                        tmp1 = A[newcand[0]][testcol] * A[leavecand[i]][col]
                        tmp2 = A[leavecand[i]][testcol] * A[newcand[0]][col]
                        if tmp1 == tmp2:
                            newcand.append(leavecand[i])
                        elif tmp1 > tmp2:
                            newcand = [leavecand[i]]
                    leavecand = newcand
                else:
                    wj = self.bascobas[j + n]
                    if wj in leavecand:
                        leavecand.remove(wj)
            if j == 0:
                z0leave = self.bascobas[0] in leavecand
            j += 1
        assert len(leavecand) == 1
        return self.whichvar[leavecand[0]], z0leave

    def negcol(self, col):
        for i in range(self.n):
            self.A[i][col] = -self.A[i][col]

    def negrow(self, row):
        for j in range(self.n + 2):
            self.A[row][j] = -self.A[row][j]

    def pivot(self, leave, enter):
        n = self.n
        A = self.A
        row = self.bascobas[leave]
        col = self.bascobas[enter] - n
        pivelt = A[row][col]
        negpiv = pivelt < 0
        if negpiv:
            pivelt = -pivelt
        for i in range(n):
            if i != row:
                nonzero = A[i][col] != 0
                for j in range(n + 2):
                    if j != col:
                        tmp1 = A[i][j] * pivelt
                        if nonzero:
                            tmp2 = A[i][col] * A[row][j]
                            if negpiv:
                                tmp1 += tmp2
                            else:
                                tmp1 -= tmp2
                        A[i][j] = tmp1 // self.determinant
                if nonzero and not negpiv:
                    A[i][col] = -A[i][col]
        A[row][col] = self.determinant
        if negpiv:
            self.negrow(row)
        self.determinant = pivelt
        self.bascobas[leave] = col + n
        self.whichvar[col + n] = leave
        self.bascobas[enter] = row
        self.whichvar[row] = enter

    def runlemke(self, *, verbose=False, lexstats=False, z0=False, silent=False):
        global filehandle
        if silent:
            filehandle = open(outfile, "w")
        n = self.n
        self.pivotcount = 1
        printout("After filltableau:")
        printout(self)
        enter = 0
        leave, z0leave = self.lexminvar(enter)
        self.negcol(n + 1)
        if verbose:
            printout("After negcol:")
            printout(self)
        while True:
            self.testtablvars()
            if z0:
                if self.bascobas[0] < n:
                    printout(
                        "step,z0=",
                        self.pivotcount,
                        self.A[self.bascobas[0]][n + 1] / self.determinant,
                    )
                else:
                    printout("step,z0=", self.pivotcount, 0.0)
            self.docupivot(leave, enter)
            self.pivot(leave, enter)
            if z0leave:
                if z0:
                    printout("step,z0=", self.pivotcount + 1, 0.0)
                break
            if verbose:
                printout(self)
            enter = self.complement(leave)
            leave, z0leave = self.lexminvar(enter)
            self.pivotcount += 1
        printout("Final tableau:")
        printout(self)
        self.createsol()
        printout(self.outsol())
        if lexstats:
            self.outstatistics()

    #######  end of class tableau


def printout(*s):
    print(*s, file=filehandle)


@click.command()
@click.argument("lcpfilename", default="lcp")
@click.option(
    "-v",
    "-verbose",
    "verbose",
    is_flag=True,
    default=False,
    help="Print intermediate tableaus.",
)
@click.option(
    "-s",
    "-silent",
    "silent",
    is_flag=True,
    default=False,
    help="Send output to <lcpfilename>.out instead of stdout.",
)
@click.option(
    "-z0",
    "z0",
    is_flag=True,
    default=False,
    help="Show value of z0 at each pivot step.",
)
def main(lcpfilename, verbose, silent, z0):
    global filehandle, outfile
    outfile = lcpfilename + ".out"
    if silent:
        filehandle = open(outfile, "w")

    printout(f"verbose={verbose} lcpfilename={lcpfilename} silent={silent} z0={z0}")
    m = lcp(lcpfilename)
    printout(m)
    printout("==================================")
    tabl = tableau(m)
    tabl.runlemke(verbose=verbose, z0=z0, silent=silent)


if __name__ == "__main__":
    main()
