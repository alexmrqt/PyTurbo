from PyTurbo import PyViterbi as viterbi

import numpy

#Quick implementation of a (7,5) convolutive code encoder
def encode75(msg):
    reg = [0,0]
    out_msg = numpy.zeros(2*len(msg), dtype=numpy.bool)

    for i in range(0, len(msg)):
        #Compute outputs
        out_msg[2*i] = msg[i]^reg[1]
        out_msg[2*i+1] = msg[i]^reg[0]^reg[1]

        #Update shift register
        reg[1] = reg[0]
        reg[0] = msg[i]

    return out_msg

#Compute branch metrics for a pair of received bits
def branch_metrics(bits_rcvd):
    ret_val = numpy.zeros(4, dtype=numpy.float32);
    cw = [[0.0,0.0], [0.0,1.0], [1.0,0.0], [1.0,1.0]] #The 4 different codewords

    bits_rcvd = numpy.array(bits_rcvd)
    cw = numpy.array(cw)

    for i in range(0, len(cw)):
        ret_val[i] = numpy.sum(numpy.abs(bits_rcvd-cw[i])**2)

    return ret_val

#Define trellis
I=2
S=4
O=4
NS = [0, 2, \
      0, 2, \
      1, 3, \
      1, 3]
OS = [0, 3, \
      3, 0, \
      1, 2, \
      2, 1]

#Create decoder instance
dec = viterbi(I, S, O, NS, OS)

#Length of the message
K_m = 100;
#Length of the coded message
K_c = 200 #Code efficiency is 1/2

#Generate message
m = numpy.random.randint(0, 2, K_m, dtype=numpy.bool)

#Encode message
c = encode75(m)

#Add noise
r = c + numpy.random.normal(0.0, 0.01, K_c)

#Compute branch metrics
bm = numpy.zeros(K_m*O, dtype=numpy.float32)
for i in range(0, K_m):
    bm[O*i:O*(i+1)] = branch_metrics([r[2*i], r[2*i+1]])

#decode message
m_hat = dec.viterbi_algorithm(K_m, -1, -1, bm);

print(numpy.array(m, dtype=numpy.int))
print('...')
print(m_hat)

##Compute BER
BER = numpy.sum(numpy.abs(m-m_hat))/K_m
print(BER)
