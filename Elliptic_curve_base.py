import math



def EEA(a,b):
    q = [0,0]
    r = [a,b]
    s = [1,0]
    t = [0,1]

    i = 1

    while r[i] > 0:
        i +=1
        q.append(r[i-2] // r[i-1])
        r.append(r[i-2] - q[i] * r[i-1])
        s.append(s[i-2] - q[i] * s[i-1])
        t.append(t[i-2] - q[i] * t[i-1])

    return r[i-1], s[i-1], t[i-1]

def invModP(x, p):
    y = x% p
    r,s,t = EEA(p,y)

    return t%p

def extended_gcd(a, b):
    if b == 0:
        return a,1,0

    d1, s1, t1 = extended_gcd(b, a % b)
    gcd = d1
    x = t1
    y = s1 - (a // b) * t1

    return gcd, x, y


class FieldNum:
    def __init__(self, value, p):
        self.value = value % p
        self.p = p

    def inversionFieldP(self, p, a):
        u, v, q = a, p, 0
        xa, xb = 1, 0
        while u != 1:
            q = v // u
            r, x = v - q * u, xb - q * xa
            v, u, xa, xb, = u, r, xa, x
        return xa % p

    def __add__(self, other):
        return FieldNum((self.value + other.value) % self.p, self.p)

    def __sub__(self, other):
        return FieldNum((self.value - other.value) % self.p, self.p)

    def __neg__(self):
        return FieldNum(-self.value, self.p)

    def __mul__(self, other):
        return FieldNum((self.value * other.value), self.p)

    def __truediv__(self, other):
        return FieldNum((self * FieldNum(invModP(other.value, self.p), self.p)).value, self.p)

    def __eq__(self, other):
        if self.value == other.value and self.p == other.p:
            return True
        return False


class FieldPoint:
    def __init__(self, curve, x, y):
        self.curve = curve
        self.x = x
        self.y = y
        self.p = curve.p

    def __add__(self, other):
        # checking identity
        if (self.x == self.curve.inf_point.x and self.y == self.curve.inf_point.y):
            return self
        elif (other.x == self.curve.inf_point.x and other.y == self.curve.inf_point.y):
            return other

        s_proj = ProjectiveFieldPoint(self.curve, FieldNum(self.x.value, self.p), FieldNum(self.y.value, self.p),
                                      FieldNum(1, self.p))
        o_proj = ProjectiveFieldPoint(other.curve, FieldNum(other.x.value, self.p), FieldNum(other.y.value, self.p),
                                      FieldNum(1, self.p))

        res_proj = (s_proj + o_proj).convertAffine()
        return res_proj

    def __sub__(self, other):
        return self + FieldPoint(self.curve, self.x, FieldNum(-self.y, self.p).value)

    def __eq__(self, other):
        if self.curve == other.curve and self.x == other.x and self.y == other.y:
            return True
        return False

    def opOrder(self, intiger):
        tmp = intiger
        ops = []
        while tmp != 1:
            if tmp % 2 == 0:
                ops.append(2)
                tmp = tmp / 2
            else:
                ops.append(1)
                tmp -= 1
        ops.reverse()
        return ops

    def simpleIntMul(self, intiger):
        temp_point = FieldPoint(self.curve, self.x, self.y)
        for i in range(intiger - 1):
            temp_point = temp_point + self
        return temp_point

    def integerMulti(self, integer):
        temp_point = FieldPoint(self.curve, self.x, self.y)
        ops = self.opOrder(integer)
        for op in ops:
            if op == 1:
                temp_point = temp_point + self
            else:
                temp_point = temp_point + temp_point

        return temp_point


class ProjectiveFieldPoint:
    def __init__(self, curve, x, y, z):
        self.curve = curve

        self.x = x
        self.y = y
        self.z = z

    def check(self):
        return self.y.value * self.y.value * self.z.value == self.x.value * self.x.value * self.x.value + self.curve.a.value * self.x.value * self.z.value * self.z.value + self.curve.b * self.z.value * self.z.value * self.z.value

    def convertAffine(self):
        return FieldPoint(self.curve, self.x / self.z, self.y / self.z)

    def __add__(self, other):
        # idea for breaking the code up into chunks comes from:
        # https://www.nayuki.io/page/elliptic-curve-point-addition-in-projective-coordinates
        if other.z.value * other.y.value == self.z.value * self.y.value:
            # same point, so need to double
            return self.double(other)
        # building points
        numeA = FieldNum(self.y.value * other.z.value, self.curve.p).value
        numeB = FieldNum(other.y.value * self.z.value, self.curve.p).value
        denomA = FieldNum(self.x.value * other.z.value, self.curve.p).value
        denomB = FieldNum(other.x.value * self.z.value, self.curve.p).value

        t = FieldNum(numeA - numeB, self.curve.p).value
        u = FieldNum(denomA - denomB, self.curve.p).value
        u2 = FieldNum(u * u, self.curve.p).value
        v = FieldNum(self.z.value * other.z.value, self.curve.p).value
        w = FieldNum(t * t * v - u2 * (denomA + denomB), self.curve.p).value
        u3 = FieldNum(u * u2, self.curve.p).value

        rx = FieldNum(u * w, self.curve.p)
        ry = FieldNum(t * (denomA * u2 - w) - numeA * u3, self.curve.p)
        rz = FieldNum(u3 * v, self.curve.p)

        return ProjectiveFieldPoint(curve=self.curve, x=rx, y=ry, z=rz)

    def __neg__(self):
        return ProjectiveFieldPoint(curve=self.curve, x=self.x, y=-self.y, z=self.z)

    def __sub__(self, other):
        return self + -other

    def __mul__(self, intiger):
        tmp = ProjectiveFieldPoint(curve=self.curve, x=self.x, y=self.y, z=self.z)

    def double(self, other):
        # idea for breaking the code up into chunks comes from:
        # https://www.nayuki.io/page/elliptic-curve-point-addition-in-projective-coordinates
        t = FieldNum(self.x.value * self.x.value * 3 + self.curve.a * self.z.value * self.z.value, self.curve.p).value
        u = FieldNum(self.y.value * self.z.value * 2, self.curve.p).value
        v = FieldNum(u * self.x.value * self.y.value * 2, self.curve.p).value
        w = FieldNum(t * t - v * 2, self.curve.p).value

        rx = FieldNum(u * w, self.curve.p)
        ry = FieldNum(t * (v - w) - u * u * self.y.value * self.y.value * 2, self.curve.p)
        rz = FieldNum(u * u * u, self.curve.p)

        return ProjectiveFieldPoint(curve=self.curve, x=rx, y=ry, z=rz)


class EllipticCurve:
    def __init__(self, a, b, mod_p):
        self.a = a
        self.b = b
        self.p = mod_p

        infNUm = FieldNum(self.p, self.p)
        infNUm.value = self.p
        self.inf_point = FieldPoint(self, infNUm, infNUm)
        self.points = [self.inf_point]

        self.gen_point = next(self.naiveGenPoints())

    def genPossibleX(self):
        for x in range(self.p):
            yield (x, FieldNum(x ** 3 + self.a * x + self.b, self.p).value)

    def genPossibleY(self):
        for y in range(self.p):
            yield (y, FieldNum(y ** 2, self.p).value)

    def naiveGenPoints(self):
        for x in self.genPossibleX():
            for y in self.genPossibleY():
                if y[1] == x[1]:
                    yield FieldPoint(self, FieldNum(x[0], self.p), FieldNum(y[0], self.p))

    def genPoints(self):
        for i in range(2, self.p):
            yield self.gen_point.integerMulti(i)

    def checkValid(self):
        return 4 * self.a ** 3 + 27 * self.b ** 2 != 0


def main():
    p = 16001
    a = 10
    b = 1
    P = (1654, 7208)
    n = 8026
    Q = (5000, 1283)

    ec = EllipticCurve(a, b, p)

    p = FieldPoint(ec, FieldNum(1654, p), FieldNum(7208, p))
    Qp = p.integerMulti(2048)
    print(Qp.x.value, Qp.y.value)





if __name__ == '__main__':
    main()
