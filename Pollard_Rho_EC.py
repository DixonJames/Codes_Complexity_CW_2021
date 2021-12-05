import logging
import math
import des
import random
import binascii
from Elliptic_curve_base import *


class Decrypt:
    def __init__(self, ciphertext, key_int):
        self.ciphertext = bytes.fromhex(ciphertext)
        self.key_int = key_int
        self.key = des.DesKey(self.deriveKey())


    def deriveKey(self):
        k56_bytes = [i for i in ("0" * (56 -len(bin(self.key_int)[2:]))) + bin(self.key_int)[2:]]
        whole = ""
        for i in range(8):
            p = k56_bytes[i * 7:(i + 1) * 7]
            tot = 0
            for c in p:
                tot += int(c)
            if tot%2 == 0:
                p.append("1")
            else:
                p.append("0")

            k8_bytes = ""
            for s in p:
                k8_bytes += s

            binary_int = int(k8_bytes, 2)
            hex_s = hex(binary_int)
            whole = whole + hex_s
            print(p)



        byte_array = bytearray(whole.encode())
        return byte_array


class PollardRho:
    def __init__(self, Q, P, p, n, curve, partitions):
        self.Q = FieldPoint(x=FieldNum(Q[0], p), y=FieldNum(Q[1], p), curve=curve)
        self.P = FieldPoint(x=FieldNum(P[0], p), y=FieldNum(P[1], p), curve=curve)
        self.prime = p
        self.n = n
        self.curve = curve

        self.part_num = partitions
        self.part_ranges = []
        self.function_coeffs = []

        self.createPartitions()

    def createPartitions(self):
        for p_n in range(self.part_num):
            a, b = self.randCoeff(), self.randCoeff()
            self.function_coeffs.append((a, b, self.P.integerMulti(a.value) + self.Q.integerMulti(b.value)))

    def randCoeff(self):
        return FieldNum(random.randint(0, self.n - 1), self.n)

    def partitionFunction(self, X):
        return int(bin(X.x.value)[-4:], 2)

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
        if self.Q == self.P.integerMulti(l):
            return True
        return False

    def coefficients(self):
        a = self.randCoeff()
        b = self.randCoeff()

        X, c, d = self.P.integerMulti(a.value) + self.Q.integerMulti(b.value), a, b
        Xp, cp, dp = self.orbitStep(X, c, d)

        exp_i = int(math.sqrt(math.pi * (self.n / 2)))
        i = 0
        while X != Xp or i == 0 or d == dp:
            i += 1
            X, c, d = self.orbitStep(X, c, d)

            Xp, cp, dp = self.orbitStep(Xp, cp, dp)
            Xp, cp, dp = self.orbitStep(Xp, cp, dp)

        if True == False:
            if self.P.integerMulti(c.value) + self.Q.integerMulti(d.value) == X and self.P.integerMulti(cp.value) + self.Q.integerMulti(
                    dp.value) == Xp and X==Xp:
                print("good coefficients")
                #print(X.x.value, X.y.value)

        return c.value, cp.value, d.value, dp.value
        # l = FieldNum((c - cp)/(dp - d), self.prime).value

    def fullPollard(self):
        c, cp, d, dp = self.coefficients()
        c, cp, d, dp = FieldNum(c, self.n), FieldNum(cp, self.n), FieldNum(d, self.n), FieldNum(dp, self.n)

        if math.gcd((dp - d).value, self.prime) == 1:
            l = ((c - cp) / (dp - d)).value
            if self.checkL(l) == True:
                return l
            else:
                return self.fullPollard()
        else:
            # here we go...
            u = FieldNum(c, self.prime) - FieldNum(cp, self.prime)
            v = FieldNum(dp, self.prime) - FieldNum(d, self.prime)

            d, x, y = extended_gcd(dp - d, self.prime)

            s = FieldNum(x, self.prime)

            w = s * u

            start = w / FieldNum(d, self.prime)
            step = FieldNum(self.n, self.prime) / FieldNum(d, self.prime)

            possible_ans = []
            for i in range(self.n - 1):
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

        # l = FieldNum((c - cp)/(dp - d), self.prime).value


def test_example():
    p = 229
    a = 10
    b = 1
    P = (5, 116)
    n = 239
    Q = (155, 166)

    ec = EllipticCurve(a, b, n - 1)

    PR = PollardRho(Q, P, p, n, ec, 16)
    c, cp, d, dp = PR.fullPollard()


def Basic_Pollard_rho():
    p = 16001
    a = 10
    b = 1
    P = (1654, 7208)
    n = 8026
    Q = (5000, 1283)

    ec = EllipticCurve(a, b, p)

    PR = PollardRho(Q, P, p, n, ec, 16)
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

    PR = PollardRho(Q, P, p, n, ec, 16)
    l = PR.fullPollard()
    print(f"l:{l}")
    return l

def decrypt():
    ciphertext = "3da46f7b6fa82f53153908bdadcc742ac38e8691e5208aa4bf6be47240c71e75180b9d1030a00810"
    l = Full_Pollard_rho()

    ws = Decrypt(ciphertext, l)




def main():
    #test_example()
    #Basic_Pollard_rho()
    #l = Full_Pollard_rho()
    decrypt()



if __name__ == '__main__':
    main()
