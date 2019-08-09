from PyTurbo import PyLogBCJR as bcjr
from PyTurbo import PyMaxLogBCJR as max_log_bcjr
from PyTurbo import PyViterbi as viterbi
from matplotlib import pyplot as plt

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

#Compute viterbi branch metrics for a sequence of received bits
def viterbi_branch_metrics(bits_rcvd, nbits_cw):
    N_cw = (2**nbits_cw)
    K = int(len(bits_rcvd)/nbits_cw)
    ret_val = numpy.zeros((K, N_cw), dtype=numpy.float32);
    cw = numpy.array([[0.0,0.0], [0.0,1.0], [1.0,0.0], [1.0,1.0]]) #The 4 different codewords

    bits_rcvd = numpy.array(bits_rcvd).reshape((K, nbits_cw))
    for k in range(0, K):
        ret_val[k][:] = numpy.sum(numpy.abs(bits_rcvd[k][:]-cw)**2, axis=1)

    return ret_val.flatten()

#Compute log-BCJR branch metrics for a sequence of received bits
def log_bcjr_branch_metrics(bits_rcvd, nbits_cw, sigma_b2):
    N_cw = (2**nbits_cw)
    K = int(len(bits_rcvd)/nbits_cw)
    ret_val = numpy.zeros((K, N_cw), dtype=numpy.float32);
    cw = numpy.array([[0.0,0.0], [0.0,1.0], [1.0,0.0], [1.0,1.0]]) #The 4 different codewords

    bits_rcvd = numpy.array(bits_rcvd).reshape((K, nbits_cw))
    for k in range(0, K):
        ret_val[k][:] = -1.0/sigma_b2 * numpy.sum(numpy.abs(bits_rcvd[k][:]-cw)**2, axis=1)

    return ret_val.flatten()

def max_log_bcjr_branch_metrics(bits_rcvd, nbits_cw, sigma_b2):
    N_cw = (2**nbits_cw)
    K = int(len(bits_rcvd)/nbits_cw)
    ret_val = numpy.zeros((K, N_cw), dtype=numpy.float32);
    cw = numpy.array([[0.0,0.0], [0.0,1.0], [1.0,0.0], [1.0,1.0]]) #The 4 different codewords

    bits_rcvd = numpy.array(bits_rcvd).reshape((K, nbits_cw))
    for k in range(0, K):
        ret_val[k][:] = -numpy.sum(numpy.abs(bits_rcvd[k][:]-cw)**2, axis=1)

    return ret_val.flatten()

#Compute bit LLR from a posteriori-probabilities
def log_bcjr_compute_llr(app, K, S):
    llr = numpy.zeros(K, dtype=numpy.float32)

    app = app.reshape((K, S, 2))
    for k in range(0, K):
        #We need copy() to make the vectors C-contiguous
        llr[k] = bcjr.max_star(app[k,:,0].copy()) - bcjr.max_star(app[k,:,1].copy())

    return llr

def max_log_bcjr_compute_llr(app, K, S):
    llr = numpy.zeros(K, dtype=numpy.float32)

    app = app.reshape((K, S, 2))
    llr = numpy.max(app[:,:,0], axis=1) - numpy.max(app[:,:,1], axis=1)
    #for k in range(0, K):
    #    llr[k] = numpy.max(app[k,:,0]) - numpy.max(app[k,:,1])

    return llr

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
R = 1/2 # Code efficiency

#Length of the message
K_m = 500000;
#Length of the coded message
K_c = int(K_m/R) #Code efficiency is 1/2

# Per-bit SNR (in dB)
EbN0dB = numpy.arange(0, 8)

# Compute noise variance from EB/N0
sigma_b2 = numpy.power(10, -EbN0dB/10)
sigma_b2 *= 1/2 # 0.5*Ps/R

#Create decoder instance
dec_vit = viterbi(I, S, O, NS, OS)
dec_log_bcjr = bcjr(I, S, O, NS, OS)
dec_max_log_bcjr = max_log_bcjr(I, S, O, NS, OS)

BER_viterbi = numpy.zeros(len(EbN0dB))
BER_log_bcjr = numpy.zeros(len(EbN0dB))
BER_max_log_bcjr = numpy.zeros(len(EbN0dB))
for i in range(0, len(EbN0dB)):
    #Generate message
    m = numpy.random.randint(0, 2, K_m, dtype=numpy.bool)

    #Encode message
    c = encode75(m)

    #Add noise
    noise = numpy.random.normal(0.0, numpy.sqrt(sigma_b2[i]), K_c) \
            + 1j*numpy.random.normal(0.0, numpy.sqrt(sigma_b2[i]), K_c)
    noise /= numpy.sqrt(2)
    r = c + noise

    ## Viterbi
    #Compute branch metrics
    bm_viterbi = viterbi_branch_metrics(r, int(1/R))

    #decode message
    m_hat_viterbi = dec_vit.viterbi_algorithm(-1, -1, bm_viterbi);

    ## Log BCJR
    #Compute branch metrics
    bm_log_bcjr = log_bcjr_branch_metrics(r, int(1/R), sigma_b2[i])

    #Compute posterior probabilities of codewords
    #A0 = numpy.log([1.0, 1e-20, 1e-20, 1e-20], dtype=numpy.float32) #Trellis begin in first state (all-0)
    A0 = numpy.log([1.0/4]*4, dtype=numpy.float32) #Do not know in which state we end
    BK = numpy.log([1.0/4]*4, dtype=numpy.float32) #Do not know in which state we end
    post_log_bcjr = dec_log_bcjr.log_bcjr_algorithm(A0, BK, bm_log_bcjr);

    #Compute bit LLR and take decisions
    llr_log_bcjr = log_bcjr_compute_llr(post_log_bcjr, K_m, S)
    m_hat_log_bcjr = (llr_log_bcjr<0)

    ## Max-Log BCJR
    #Compute branch metrics
    bm_max_log_bcjr = max_log_bcjr_branch_metrics(r, int(1/R), sigma_b2[i])

    #Compute posterior probabilities of codewords
    post_max_log_bcjr = dec_max_log_bcjr.log_bcjr_algorithm(A0, BK, bm_max_log_bcjr);

    #Compute bit LLR and take decisions
    llr_max_log_bcjr = max_log_bcjr_compute_llr(post_max_log_bcjr, K_m, S)
    m_hat_max_log_bcjr = (llr_max_log_bcjr<0)

    ##Compute BER
    BER_viterbi[i] = numpy.mean(numpy.abs(m!=m_hat_viterbi))
    print('BER For viterbi at Eb/N0 = ' + str(EbN0dB[i]) + 'dB: ' + str(BER_viterbi[i]))
    BER_log_bcjr[i] = numpy.mean(numpy.abs(m!=m_hat_log_bcjr))
    print('BER For log_bcjr at Eb/N0 = ' + str(EbN0dB[i]) + 'dB: ' + str(BER_log_bcjr[i]))
    BER_max_log_bcjr[i] = numpy.mean(numpy.abs(m!=m_hat_max_log_bcjr))
    print('BER For max_log_bcjr at Eb/N0 = ' + str(EbN0dB[i]) + 'dB: ' + str(BER_max_log_bcjr[i]))
    print('')

plt.semilogy(EbN0dB, BER_viterbi, label="Viterbi")
plt.semilogy(EbN0dB, BER_log_bcjr, label="Log-BCJR")
plt.semilogy(EbN0dB, BER_max_log_bcjr, label="Max-Log-BCJR")
plt.legend()
plt.xlabel("E_b/N_0 (dB)")
plt.ylabel("BER")
plt.grid(True)
plt.show()
