from lab2code import Question_1b
import numpy as np
import matplotlib.pyplot as plt
import math as m

def Question_3b(Filename):
    grad, array = Question_1b(Filename)                                 # read data from the text file
    X, Y = np.meshgrid(np.arange(0, 1201, 10), np.arange(0, 1201, 10))
    plt.figure()
    plt.quiver(X, Y, grad[0], grad[1])
    plt.title("The gradient map of Big Island of Hawaii")
    plt.xlabel("X data sample points")
    plt.ylabel("Y data sample points")
    plt.savefig("gradient_map.png")
    plt.show()

def gradient_real(data):
    transform = np.fft.fft(data)
    k = np.fft.fftfreq(1201)
    diff_coeff = 1j * 2 * m.pi * k
    ft_diff = transform * diff_coeff
    grad = np.fft.ifft(ft_diff)
    return grad.real

def Question_3c(Filename):
    grad, array = Question_1b(Filename)         # read data from the text file
    dx = gradient_real(array)
    array = np.transpose(array)
    dy = gradient_real(array)
    plt.imshow(dx,  vmin=-0.1, vmax=0.075)      # choose max and min for the height
    plt.title("Gradient map generate with FFT.")
    plt.xlabel("x data points")
    plt.ylabel("y data points")
    plt.savefig("FFT generate.png")
    plt.show()






if __name__ == "__main__":
    Question_3b("N19W156.hgt")
    Question_3c("N19W156.hgt")
