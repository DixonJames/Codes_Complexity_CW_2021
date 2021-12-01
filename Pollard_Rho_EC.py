import math
import random
from Elliptic_curve_base import *

class PollardRho:
    def __init__(self, Q, P, p, n, curve, partitions):
        self.Q = FieldPoint(x=FieldNum(Q[0], p), y=FieldNum(Q[1],p), curve=curve)
        self.P = FieldPoint(x=FieldNum(P[0], p), y=FieldNum(P[1], p), curve=curve)
        self.prime = p
        self.n = n
        self.curve = curve

        self.part_num = partitions
        self.part_ranges = []
        self.function_coeffs = []

        self.createPartitions()

    def createPartitions(self):
        part_size = math.floor(self.prime/self.part_num)
        diff = self.prime - part_size * self.part_num  # just add to last part

        c = 0
        for p_n in range(self.part_num):
            self.part_ranges.append((p_n * part_size, (p_n + 1) * part_size -1))
            self.function_coeffs.append((self.randCoeff(), self.randCoeff()))

        self.part_ranges[-1] = (self.part_ranges[-1][0], self.part_ranges[-1][1] + diff + 1)

    def randCoeff(self):
        return random.randint(0, self.n - 1)

    def partitionFunction(self, X):
        x = X.x

        c = 1
        i = len(self.part_ranges)//2
        found = False
        while not found:
            c += 1
            c_range = self.part_ranges[i]
            if x.value < c_range[0]:
                i = i - math.ceil(len(self.part_ranges) * (1/2) ** c)
            elif x.value > c_range[1]:
                i = i - math.ceil(len(self.part_ranges) * (1/2) ** c)
            else:
                return i

    def orbitStep(self, X,c,d):
        S_num = self.partitionFunction(X)
        a = (self.function_coeffs[S_num][0]) #% (self.n - 1)
        b = (self.function_coeffs[S_num][1]) #% (self.n - 1)
        new_a = self.function_coeffs[S_num][0] + a
        new_b = self.function_coeffs[S_num][1] + b
        return X + self.P.integerMulti(a) + self.Q.integerMulti(b), new_a, new_b

    def BasicPollard(self):
        a = self.randCoeff()
        b = self.randCoeff()
        X, c, d = self.P.integerMulti(a) + self.Q.integerMulti(b), a, b
        Xp, cp, dp = self.orbitStep(X, a, b)
        cp += c
        dp += d

        i = 0
        while X != Xp or (cp == c or dp == d):
            i+= 1
            X, c, d = self.orbitStep(X, c, d)
            """c = (c + ct) #% (self.n - 1)
            d = (d + dt) #% (self.n - 1)"""

            Xp, cp, dp = self.orbitStep(Xp, cp, dp)
            """cp = (cp + cpt) #% (self.n - 1)
            dp = (dp + dpt) #% (self.n - 1)"""

            Xp, cp, dp = self.orbitStep(Xp, cp, dp)
            """cp = (cp + cpt) #% (self.n - 1)
            dp = (dp + dpt) #% (self.n - 1)"""

            #print(c, cp, d, dp)

        return c, cp, d, dp
        #l = FieldNum((c - cp)/(dp - d), self.prime).value

    def checkL(self, l):
        if self.Q == self.P.integerMulti(l):
            return True
        return False

    def fullPollard(self, a = 5, b = 10):
        X, c, d = self.P.integerMulti(a) + self.Q.integerMulti(b), a, b
        Xp, cp, dp = self.orbitStep(X)

        while X != Xp:
            X, ct, dt = self.orbitStep(X)
            c += ct
            d += dt

            Xp, cpt, dpt = self.orbitStep(Xp)
            cp += cpt
            dp += dpt

            Xp, cpt, dpt = self.orbitStep(Xp)
            cp += cpt
            dp += dpt

            #print(c, cp, d, dp)

        if math.gcd(dp - d, self.prime) == 1:
            l = ((FieldNum(c, self.prime) - FieldNum(cp, self.prime))/(FieldNum(dp, self.prime) - FieldNum(d, self.prime))).value
            if self.checkL(l):
                return l
        else:
            #here we go...
            u = FieldNum(c, self.prime) - FieldNum(cp, self.prime)
            v = FieldNum(dp, self.prime) - FieldNum(d, self.prime)

            d, x, y = extended_gcd(dp - d, self.prime)

            s = FieldNum(x, self.prime)

            w = s * u

            start = w / FieldNum(d, self.prime)
            step = FieldNum(self.n, self.prime) / FieldNum(d, self.prime)

            possible_ans = []
            for i in range(self.n-1):
                possible_ans.append((start + FieldNum(step.value * i, self.prime)).value)

            for l in possible_ans:
                if self.checkL(l):
                    return l

        """
        if d == 1: return ((k[2] - j[2]) * inverse(j[1] - k[1], p - 1)) % (p - 1)
        m, l = 0, ((k[2] - j[2]) * inverse(j[1] - k[1], (p - 1) / d)) % ((p - 1) / d)
        while m <= d:
            print
            m, l
            if pow(g, l, p) == t: return l
            m, l = m + 1, (l + ((p - 1) / d)) % (p - 1)
        return False
        """


        #l = FieldNum((c - cp)/(dp - d), self.prime).value

def Basic_Pollard_rho():
    p = 16001
    a = 10
    b = 1
    P = (1654, 7208)
    n = 8026
    Q = (5000, 1283)

    ec = EllipticCurve(a, b, n-1)

    PR = PollardRho(Q, P, p, n, ec, 50)
    c, cp, d, dp = PR.BasicPollard()
    print(f"c:{c}, c':{cp}, d:{d}, d':{dp}")



def main():
    p = 16001
    a = 10
    b = 1
    P = (1654, 7208)
    n = 8026
    Q = (5000, 1283)

    ec = EllipticCurve(a, b, p)

    PR = PollardRho(Q, P, p, n, ec, 3)
    l = PR.BasicPollard()
    print(f"l:{l}")



if __name__ == '__main__':
    main()
