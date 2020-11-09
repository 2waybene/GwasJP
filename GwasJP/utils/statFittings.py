import math
import numpy


def gaussian(x, mu, sig):
    '''
    Gaussian function.
    '''
    return (1 / math.sqrt(2 * math.pi * sig ** 2)) * \
        numpy.exp(-(x - mu) ** 2 / (2 * sig ** 2))


def logistic(x, x0, L, M, k):
    '''
    Logistic function.
    '''
    return M + (L / (1 + numpy.exp(-k * (x - x0))))

def entropy (num1, num2):
    percnt1 = num1/(num1+num2)
    percnt2 = num2/(num1+num2)
    return (-(percnt1*math.log2(percnt1) + percnt2*math.log2(percnt2)))



if __name__ == "__main__":
    print(entropy(5,9))
    print(entropy(3,2))
    print(entropy(2,3))
