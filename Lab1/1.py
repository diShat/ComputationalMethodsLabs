# Calc Methods Lab 1 (Solving Nonlinear Equations) by Diana Shatkovska FB-92, Variant 3

accuracy = 0.00001
p_coef = (1, -3, 1, -2, -4)
intervals = [(-2.41, -0.57), (0.66, 5)]


def f(x, coef=p_coef):
    res = 0
    for i, k in enumerate(coef[::-1]):
        res += k * x**i

    return res


def f_derivative(x, coef=p_coef):
    res = 0
    for i, k in enumerate(coef[::-1]):
        res += (k*i) * x**(i-1)

    return res


def bisection(interval, e=accuracy):
    print('---Bisection method---')
    left, right = interval
    print(f"0:\t{left},{right}")
    i = 0
    while abs(right-left) > e:
        bisector = (right + left)/2

        if f(left)*f(bisector) <= 0:
            right = bisector
        else:
            left = bisector
        i += 1
        print(f"{i}:\tx = {(right + left)/2}, e = {abs(right-left):.6f}")
    print(f"Result:\tx = {(right + left)/2}")
    return (right + left)/2


# bisection(intervals[0])
# bisection(intervals[1])


def secant(interval, e=accuracy):
    print('---Secant method---')
    span = interval    # c, c_prev
    print(f"0:\t{span}")
    i = 0

    while abs(span[1]-span[0]) > e or abs(f(span[0])) > e:
        span = (span[0] * f(span[1]) - span[1] * f(span[0])) / (f(span[1])-f(span[0])), span[0]
        i += 1
        print(f"{i}:\tx = {span[0]}, e = {abs(span[1]-span[0]):.6f}")
    print(f"Result:\tx = {span[0]}")
    return span[0]


# secant(intervals[0])
# secant(intervals[1])
# secant((2.125,5))


def newtons(interval, e=accuracy):
    print('---Newton method---')
    span = interval     # x, x0
    print(f"0:\t{span}")
    i = 0

    while abs(span[1]-span[0]) > e and abs(f(max(span))) > e:
        span = span[1] - f(span[1])/f_derivative(span[1]), span[0]

        i += 1
        print(f"{i}:\tx = {span[0]}, |b-a| = {abs(span[1]-span[0]):.6f}, |f(a)| = {abs(f(max(span))):.6f}")
    print(f"Result:\tx = {span[0]}")
    return span[0]


# newtons(intervals[0])
# newtons(intervals[1])
# newtons((2.125,5))