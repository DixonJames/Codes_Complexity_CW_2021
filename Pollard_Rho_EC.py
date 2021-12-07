import logging
import math
import des
import random
import binascii
from Elliptic_curve_base import *


class Decrypt:
    def __init__(self, ciphertext, key_int):
        self.c = ciphertext
        self.ciphertext = bytes.fromhex(ciphertext)
        self.key_int = key_int
        self.key = des.DesKey(self.deriveKey())

        self.plaintext = self.decrypt().decode('uft-8')

    def deriveKey(self):
        b_rep = bin(self.key_int)[2:]
        k56_bytes = "0" * (56 - len(b_rep)) + b_rep
        whole = ""
        for i in range(8):
            p = k56_bytes[i * 7:(i + 1) * 7]
            tot = 0
            for c in p:
                tot += int(c)
            if tot % 2 == 0:
                p += ("1")
            else:
                p += ("0")

            whole = whole + p

        bytes_8 = int(whole, 2).to_bytes((len(whole) + 7) // 8, byteorder='big')
        return bytes_8

    def decrypt(self):
        return self.key.decrypt(self.ciphertext, padding=True)


class PollardRho:
    def __init__(self, Qa, Qb, P, p, n, curve, partitions):
        self.Qa = FieldPoint(x=FieldNum(Qa[0], p), y=FieldNum(Qa[1], p), curve=curve)
        self.Qb = FieldPoint(x=FieldNum(Qb[0], p), y=FieldNum(Qb[1], p), curve=curve)
        self.P = FieldPoint(x=FieldNum(P[0], p), y=FieldNum(P[1], p), curve=curve)
        self.prime = p
        self.n = n
        self.curve = curve

        self.part_num = partitions
        self.part_ranges = []
        self.function_coeffs = []

        self.createPartitions()

        self.l = self.fullPollard()

        self.secret = self.sharedSecret()

    def sharedSecret(self):
        return self.Qb.integerMulti(self.l).x.value

    def createPartitions(self):
        self.function_coeffs = []
        for p_n in range(self.part_num):
            a, b = self.randCoeff(), self.randCoeff()
            self.function_coeffs.append((a, b, self.P.integerMulti(a.value) + self.Qa.integerMulti(b.value)))

    def randCoeff(self):
        return FieldNum(random.randint(0, self.n - 1), self.n)

    def partitionFunction(self, X):
        s = bin(X.x.value)[-4:]
        if X.x.value < 8:
            s = "0" * (4 - len(bin(X.x.value)[2:])) + bin(X.x.value)[2:]
        try:
            a = int(s, 2)
        except:
            print("d")
        return a

    def orbitStep(self, X, c, p):
        S_num = self.partitionFunction(X)
        a = self.function_coeffs[S_num][0]
        b = self.function_coeffs[S_num][1]
        R = self.function_coeffs[S_num][2]

        new_a = (c + a)
        new_b = (p + b)
        new_x = X + R
        return new_x, new_a, new_b

    def checkL(self, l):
        if self.Qa == self.P.integerMulti(l):
            return True
        return False

    def coefficients(self):
        a = self.randCoeff()
        b = self.randCoeff()

        X, c, d = self.P.integerMulti(a.value) + self.Qa.integerMulti(b.value), a, b
        Xp, cp, dp = self.orbitStep(X, c, d)

        exp_i = int(math.sqrt(math.pi * (self.n / 2)))
        i = 0
        while X != Xp or i == 0 or d == dp:
            i += 1
            X, c, d = self.orbitStep(X, c, d)

            Xp, cp, dp = self.orbitStep(Xp, cp, dp)
            Xp, cp, dp = self.orbitStep(Xp, cp, dp)

            if i % 10000 == 0:
                print(f"~{100 * i / exp_i}% done")

        if True == False:
            if self.P.integerMulti(c.value) + self.Qa.integerMulti(d.value) == X and self.P.integerMulti(
                    cp.value) + self.Qa.integerMulti(
                    dp.value) == Xp and X == Xp:
                print("good coefficients")
                # print(X.x.value, X.y.value)

        return c.value, cp.value, d.value, dp.value
        # l = FieldNum((c - cp)/(dp - d), self.prime).value

    def gcd_m_a(self, c, cp, d, dp):
        cpp = c - cp
        dpp = dp - d

        g, x, y = extended_gcd(dpp.value, self.n)

        new_mod = int(self.n / g)
        c_base = FieldNum(cpp.value / g, self.n)
        b_base = FieldNum(dpp.value / g, self.n)

        base_point = (c_base / b_base).value

        step = FieldNum(new_mod, self.n).value

        possible_ans = []
        for i in range(g):
            t = base_point
            for j in range(i):
                t = t + step
            possible_ans.append(t)

        return possible_ans

    def gcd_m_b(self, c, cp, d, dp):
        cpp = c - cp
        dpp = dp - d

        g, x, y = extended_gcd(dpp.value, self.n)

        v = FieldNum((x * dpp.value) % self.n, self.n).value
        w = FieldNum((x * cpp.value) % self.n, self.n).value

        possible_ans = []
        for k in range(self.n):
            t = (w / v + k * (self.n / v))%(self.n/g)
            possible_ans.append(t)
        print(set(list(set(possible_ans))))
        return list(set(possible_ans))

    def fullPollard(self):
        c, cp, d, dp = self.coefficients()
        c, cp, d, dp = FieldNum(c, self.n), FieldNum(cp, self.n), FieldNum(d, self.n), FieldNum(dp, self.n)

        if math.gcd((dp - d).value, self.n) == 1:
            l = ((c - cp) / (dp - d)).value
            if self.checkL(l) == True:
                return l
                """elif self.checkL(l / 2) == True:
                    return l / 2"""
            else:
                self.createPartitions()
                return self.fullPollard()

        else:
            l_a = self.gcd_m_a(c, cp, d, dp)
            l_b= self.gcd_m_b(c, cp, d, dp)
            # here we go...


            for l in l_a:
                if self.checkL(l):
                    return l
            for l in l_b:
                if self.checkL(l):
                    return l

            self.createPartitions()
            return self.fullPollard()

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

        # l = FieldNum((c - cp)/(dp - d), self.prime).value


def test_example():
    p = 229
    a = 10
    b = 1
    P = (5, 116)
    n = 239
    Q = (155, 166)

    ec = EllipticCurve(a, b, n - 1)

    PR = PollardRho(Q, Q, P, p, n, ec, 16)
    c, cp, d, dp = PR.fullPollard()


def Basic_Pollard_rho():
    p = 16001
    a = 10
    b = 1
    P = (1654, 7208)
    n = 8026
    Q = (5000, 1283)

    ec = EllipticCurve(a, b, p)

    PR = PollardRho(Q, Q, P, p, n, ec, 16)
    for i in range(10):
        c, cp, d, dp = PR.coefficients()
        print(f"c:{c}, c':{cp}, d:{d}, d':{dp}")


def Full_Pollard_rho():
    p = 16001
    a = 10
    b = 1
    P = (1654, 7208)
    n = 8026
    Q = (5000, 1283)

    ec = EllipticCurve(a, b, p)

    PR = PollardRho(Q, Q, P, p, n, ec, 16)
    l = PR.l
    print(f"l:{l}")
    return l


def decrypt():
    ciphertext = "3da46f7b6fa82f53153908bdadcc742ac38e8691e5208aa4bf6be47240c71e75180b9d1030a00810"

    p = 20376993552394903
    a = 10
    b = 1
    P = (1983, 6761152449250519)
    n = 1852453970120513
    QA = (18586784116581871, 12161036958498472)
    QB = (18432261261031243, 11140924411855488)

    ec = EllipticCurve(a, b, p)

    """PR = PollardRho(QB, QA, P, p, n, ec, 16)
    secret = PR.secret
    print(PR.l, secret)"""

    l_atob = 1682779984167835
    secret = 6714934996831608

    ws = Decrypt(ciphertext, secret)
    print(ws.plaintext, ws.key)


def main():
    # test_example()
    # Basic_Pollard_rho()
    for _ in range(1000):
        Full_Pollard_rho()
    #
    # decrypt()


if __name__ == '__main__':
    main()
