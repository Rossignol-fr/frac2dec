# frac2dec

This python file allow the use of exact fractional computation and rational approximation with asked precision

examples:

frac2dec('11/17) produce the periodical scripture of this rational number: '0.[6470588235294117]'

dec2frac('-12.3456[789]') produce the fractional scripture of the decimal number with infinite numbers :-12.3456789789789789789789... :
'-41111111/3330000'

others stranges number scriptures like 9.9[9] is interpreted automatically as the number 10 or 0.0000[9] by the number 0.0001 thanks to an improvement in explose_dec() function

others python functions  allow to obtain the continued fraction for a decimal number and get the differents residus:
fracont('12.829[3]') produce '[12;1, 4, 1, 6, 9]'
reduite('[12;1, 4, 1, 6, 9]') produce ['12', '13', '64/5', '77/6', '526/41', '4811/375']

other functions exist like racinek(a, n, k ) which allows to compute n decimales of the kth rooth of a given fraction 'a'

