# Outils de calculs avec les nombres réels:
# écriture des rationnels, approximations d'une fraction
# approximation des réels, fractions continues

def explose_frac(r):
    """renvoie le numérateur et le dénominateur d'une fraction écrite sous
    forme de chaîne de caractères ainsi qu'un flag d'approximation"""
    assert isinstance(r, str) and '.' not in r,\
    "explose_frac attends une chaîne de caractères sans . "
    apx = False
    if r[0]=='~': #La fraction est issue d'un calcul approché
        r = r[1:] #on retire le flag
        apx = True
    if '/' in r:
        pd = r.index('/') #position diviseur
        return int(r[:pd]), int(r[pd+1:]), apx
    else:
        return int(r), 1, apx

def explose_dec(r):
    """Transforme un décimal écrit sous forme de chaîne de caractères
    en liste sous la forme: [sgn,pE,pF,p, apx] où sgn est le signe (1 ou -1)
    pE la partie Entière, pF la partie décimale fixe, p la partie décimale
    périodique et apx un booléen indiquant si c'est une valeur approchée.
    Simplifie les cas farfelus p=9 ou p=0"""
    sgn = 1 # pour le signe
    apx = False
    if r[0] == '~' : #Le nombre décimal est une valeur approchée
        r = r[1:] #on retire le flag
        apx = True
    if r[0] == '-':
        sgn = -1
        r = r[1:] #on retire le signe
    if '.' not in r: # c'est un entier
        return sgn, r, '', '', apx
    else:
        pv = r.index('.') # donne la position du séparateur décimal
        pE = r[:pv]
    if '[' not in r:
        return sgn, pE, r[pv+1:], '', apx
    else: #Le nombre est en écriture décimale périodique
        pf = r.index('[') # pour la partie décimale périodique
        pF = r[pv+1:pf]
        pr = r[pf+1:-1]
        if pr == '0': #Traite le cas inutile x.x[0]
            pr = ''
        elif pr == '9': # Traite les cas farfelus x.x[9]
            pr = ''
            sn = pE + pF
            lsn = len(sn)
            sn1 = str(int(sn)+1)
            bz = '0'*(lsn - len(sn1))
            sn = bz + sn1
            if len(sn) == lsn + 1:
                if len(bz)>0:
                    sn=sn[1:]
                pv+=1
            lf = pf-pv
            while lf>0 and sn[-1]=='0':
                sn = sn[:-1]
                lf-=1
            pE = sn[:pv]
            pF = sn[pv:]
        return sgn, pE, pF, pr, apx

def frac2dec(r, n=-1):
    """Donne l'écriture décimale de la fraction p/q écrite
    sous la forme d'une chaîne de caractères. Si n est renseigné donne
    n chiffres après la virgule (par arrondi)
    Exemples:
    frac2dec('1/17') renvoie '0.[0588235294117647]'
    frac2dec('1/30') renvoie '0.0[3]'
    frac2dec('1/256') renvoie '0.00390625'
    frac2dec('17/12',5) renvoie '~1.41667'
    frac2dec('12/17',5) renvoie '~0.70588' """
    p, q, apx = explose_frac(r)
    sg, sgn = '', 1
    if p < 0:
        sg, sgn = '-', -1
        p = abs(p)
    r = p%q
    Lr = [r] #Liste des restes
    if r == 0: # C'est un entier !
        return frc((sgn*p//q ,1 ,apx))
    fs = sg + str(p//q) + '.' #formatage de la valeur décimale
    ls = len(fs)
    if n == 0: #on ne veut que la partie entière
        return chr(126) +fs[:-1]
    elif n == -1: #on veut tout le développement décimal
        j = n
        flg = False
        ic = 0
    else: # On ne veut que n chiffres après la virgule
        j = -n - 1 #On va chercher une décimale supplémentaire pour l'arrondi
        flg = True
        ic = 1
    while j<0:
        j += ic
        fs = fs + str(10*r//q)
        r = 10*r%q
        if n==-1 and r in Lr:
            idx = Lr.index(r)
            fs = fs[:ls+idx] + '[' + fs[ls+idx:] + ']'
            flg = False
            break
        Lr.append(r)
        if r == 0 and j<0: # c'est un décimal !
            flg = False
            break
    if flg:
    #si on est là c'est qu'on attend une valeur approchée
        fs = arrondi(fs,n)
    if apx and '~' not in fs:
        return '~'+fs
    else:
        return fs


def simplifieD(r):
    """Simplifie si possible l'écriture du décimal écrit sous forme de
    chaîne de caractères. Par exemple: 2.26000 ou 2."""
    if '.' in r:
        idc = r.index('.')
        d = r[idc+1:]
        if len(d)>0:
            while len(d)>0 and d[-1]=='0':
                d = d[:-1]
            if len(d)==0:
                return r[:idc]
            else:
                return r[:idc+1]+d
        else:
            return r[:idc]
    else:
        return r

def simpfrac(a,b,apx):
    """Simplifie la fraction a/b en déterminant le pgcd,
    préserve le signe et l'ordre en renvoyant (num, den)"""
    swap = False
    if a*b < 0: # sauvegarde du signe
        sgn = -1
    else:
        sgn = 1
    a, b = abs(a), abs(b)
    if a<b:
        a, b = b, a
        swap = True
    d = pgcd(a,b)
    a, b = a//d, b//d
    if swap:
        return sgn*b, a, apx
    else:
        return sgn*a, b, apx

def pgcd(a, b):
    """Algorithme d'Euclide sur les entiers naturels a et b
    où a>b et renvoie le dernier reste non nul"""
    while a%b != 0:
        a, b = b, a%b
    return b

def dec2frac(r):
    """Converti un décimal écrit sous forme de chaîne de caractères
    en fraction (string)"""
    sgn, pE, pF, p, apx = explose_dec(r)
    t = len(p) #longueur de la période
    if t > 0:
        s = len(pF) # nombre de décimales non périodiques
        dd = pE + pF + p
        d = pE + pF
        f = simpfrac(int(dd)-int(d), sgn*(10**t-1)*10**s, apx)
    else:
        nd = len(pF)
        f = simpfrac(int(pE)*10**nd+int(pF), sgn*10**nd, apx)
    return frc(f)

def frc(f):
    """transforme f = (p, q, apx) en fraction (str)"""
    ap=''
    if f[2]:
        ap+='~'
    if f[1] == 1:
        return ap + str(f[0])
    else:
        return ap + str(f[0]) + '/' + str(f[1])

def prodf(r1, r2, simp = True):
    """renvoie le produit de deux fractions écrites avec le symbole /"""
    p1, q1, apx1 = explose_frac(r1)
    p2, q2, apx2 = explose_frac(r2)
    apx = apx1 or apx2
    if simp:
        return frc(simpfrac(p1*p2, q1*q2, apx))
    else:
        return frc((p1*p2, q1*q2, apx))

def powf(r,n):
    """renvoie la fraction élevée à la puissance n"""
    if n==1:
        return r
    elif n%2==0:
        return powf(prodf(r,r,False), n//2)
    else:
        return prodf(r,powf(prodf(r,r,False), (n-1)//2))

def addf(r1, r2, simp = True):
    """renvoie la somme de deux fractions écrites avec le symbole /"""
    p1, q1, apx1 = explose_frac(r1)
    p2, q2, apx2 = explose_frac(r2)
    apx = apx1 or apx2
    if p1*q2 + q1*p2 == 0:
        return '0'
    else:
        if simp:
            return frc(simpfrac(p1*q2 + q1*p2, q1*q2, apx))
        else:
            return frc((p1*q2 + q1*p2, q1*q2, apx))

def invf(r):
    p, q, apx = explose_frac(r)
    return frc((q, p, apx))

def negf(r):
    p, q, apx = explose_frac(r)
    if p*q<0:
        return frc((abs(p), abs(q), apx))
    else:
        return frc((-abs(p), abs(q), apx))

def racEk(m, k = 2):
    """Détermine la racine kième entière de l'entier m par dichotomie"""
    t = len(str(m))-1
    a = 10**(t//k)
    b= 10**(t//k+1)-1
    while b-a>1:
        d = (a+b)//2
        if d**k == m:
            return d
        elif d**k < m:
            a = d
        else:
            b = d
    if b**k-m < m-a**k:
        return b
    else:
        return a

def racinek(a, n, k = 2):
    """Renvoie le décimal approchant la racine kième de a avec n décimales où
    a est une fraction (str)
    Exemple du calcul de la racine cubique de 1/2:
    racinek('1/2',20,3) renvoie '~0.79370052598409973738' """
    p, q, apx = explose_frac(a)
    s = str(racEk(p*10**(n*k + k), k)) + '/' + str(racEk(q*10**(n*k + k), k))
    v = frac2dec(s, n)
    return v

def quotient_10n(s, n):
    """Détermine si |s| > 10^-n où s est une fraction (string)
    Exemple quotient_10n('1/512',3) renvoie True car 1/512>1/1000"""
    p, q, apx = explose_frac(s)
    if p == 0:
        return False
    if abs(q//p) < 10**n:
        return True
    else:
        return False

def arrondi(sd, n):
    """Renvoie l'arrondi à n décimales du nombre décimal (string) en
    prenant en compte l'écriture décimale périodique.
    Exemple: arrondi('-1.23[4567890]',4) renvoie: '-1.2346' """
    sgn, pE, pF, p, apx = explose_dec(sd)
    sg, ap = '', ''
    if sgn == -1:
        sg = '-'
    if apx:
        ap='~'
    if len(p)>0:
        while len(pF) < n+1 :
            pF+=p
    elif len(pF) <= n: #Dans ce cas le décimal est déjà formaté
        return simplifieD(ap + sg + pE +'.'+ pF)
    sn = pE + pF[:n + 1]
    pv = len(pE) # position de la virgule
    nc = len(sn)
    nci = len(str(int(sn)))
    nz = nc - nci # cas où on a 0.xx
    bz = '0'*nz # pour le bourrage de zéros
    sn = str(int(sn)+5)
    if len(sn) > nci:
        if len(bz)>0:
            bz=bz[1:]
        else:
            pv+=1
    ap = '~'
    sn = bz + sn
    return simplifieD( ap + sg + sn[:pv] + '.' + sn[pv:-1])

def approxE(n):
    """Approximation de e (exp(1)) par la série somme(1/k!,k ) où le dernier
    élément ajouté est le premier inférieur ou égal à 10^-n"""
    p, fac, f = 1, 1, '1'
    while fac < 10**n:
        f = addf(f,frc((1, fac, False)))
        p += 1
        fac = fac*p
    return f

def approxPi2s6(n, simp = True):
    """serie des inverses des carrés (converge très très lentement)
    n représente le dernier résidu inférieur à 10^-n"""
    p, f = 1, '1'
    t, T = 0, True
    j = 5
    r = 0
    while p*p <= 10**(n+j):
        p += 1
        f = addf(f,frc((1, p*p, False)), simp)
        if p*p > 10**(n+1):
            t += 1
            if p*p > 10**(n+4):
                r += 1
    return f, t, r

def lnp(n,p):
    """renvoie une fraction comme approximation du logarithme népérien
     de l'entier p grâce à une série à n termes"""
    f = frc((p-1, p+1,True))
    for k in range(1, n):
        f = addf(f, frc(((p-1)**(2*k+1), (2*k+1)*(p+1)**(2*k+1),True)), True)
    return prodf('2', f)

def fracont(d, n = -1):
    """Renvoie la fraction continue du nombre d écrit en chaîne de caractères
    d peut être soit un décimal soit une fraction. Si n!=1, donne n résidus
    Exemples:
    fracont('103/101') renvoie: '[1; 50, 2]'
    fracont('-34/21') renvoie '-[1; 1, 1, 1, 1, 1, 2]'"""

    if '.' in d:
        p, q, apx = explose_frac(dec2frac(d))
    else:
        p, q, apx = explose_frac(d)
    sg = ''
    if p*q < 0:
        sg = '-'
    p, q = abs(p), abs(q)
    pE = p//q
    s = sg + '[' + str(pE) + '; '
    p -= pE*q
    j = -n
    ic = 1
    if n == -1:
        ic = 0
        j = -1
    while j < 0 and p != 0:
        p, q = q, p
        r = p // q
        s = s+  str(r)+ ', '
        p -= r*q
        j += ic
    s = s[:-2]
    return s + ']'

def reduite(fc,n=-1):
    """Renvoie la réduite de rang n ou une liste contenant
    toutes les réduites si n=-1"""
    sg = ''
    Lr = []
    if fc[0] == '-':
        sg = '-'
        fc = fc[1:]
    idx = fc.find(';')
    pE = fc[1:idx]
    Lr.append(sg+pE)
    if idx == len(fc)-2 or n == 0:
        return Lr[0]
    fc = fc[idx+2:-1]
    La = fc.split(', ')
    Li = [int(k) for k in La]
    p0 ,p1, q0, q1 = 0, 1, 1, Li[0]
    Lr.append(sg+addf(pE,str(p1)+'/'+str(q1)))
    if n == 1 or len(Li) == 1:
        return Lr[1]
    else:
        t = n
        if n> len(Li) or n==-1:
            t = len(Li)
        for k in range(1,t):
            p = Li[k]*p1 + p0
            q = Li[k]*q1 + q0
            p0, q0 = p1, q1
            p1, q1 = p, q
            Lr.append(sg+addf(pE,str(p1)+'/'+str(q1)))
        if n == -1:
            return Lr
        else:
            return Lr[-1]

#--------------- Tests --------------
# ils mettent en avant des cas limites ou des cas standards d'erreur
#
# Test des fcts arrondi et explose_dec

dec1 = '1.23[45]'
dec2 = '0.125'
dec3 = '~-0.00095'
dec4 = '-0.999'
dec5 = '1.[9]'
dec6 = '-99.[9]'
dec7 = '0.00009[9]'
dec8 = '12.23[0]'
dec9 = '-0.00042[9]'

assert explose_dec(dec1)==(1, '1', '23', '45', False), 'Pb avec explose_dec'
assert explose_dec(dec3)==(-1, '0', '00095', '', True), 'Pb avec explose_dec'
assert explose_dec(dec6)==(-1, '100', '', '', False), 'Pb avec explose_dec'
assert explose_dec(dec8)==(1, '12', '23', '', False), 'Pb avec explose_dec'
assert explose_dec(dec7)==(1, '0', '0001', '', False), 'Pb avec explose_dec'

assert arrondi(dec1,3)=='~1.235', "Problème d'arrondi"
assert arrondi(dec1,0)=='~1', "Arrondi entier erroné"
assert arrondi(dec2,3)==dec2, "Arrondi erroné"
assert arrondi(dec3,3)=='~-0.001', "Erreur d'arrondi !"
assert arrondi(dec3,4)=='~-0.001', "Erreur d'arrondi !"
assert arrondi(dec4,2)=='~-1', "Arrondi erroné"
assert arrondi(dec5,0)=='2', "Arrondi erroné"
assert arrondi(dec6,0)=='-100', "Arrondi erroné"


# Test de la fonction frac2dec
fr1 = '1/8'
fr2 = '-12/4'
fr3 = '1/17'
fr4 = '77708431/2640858'

assert frac2dec(fr1,3)=='0.125' , 'problème dans frac2dec'
assert frac2dec(fr2,1)=='-3' , 'problème dans frac2dec'
assert frac2dec(fr3,8)=='~0.05882353' , 'problème dans frac2dec'
assert frac2dec(fr3)=='0.[0588235294117647]','problème dans frac2dec'

#n=200
#print("ln(2)",frac2dec(ln2(n),n))

#print("sqrt(2)",frac2dec(extracine('2',n)[0],n))
print(reduite(fracont(prodf(lnp(100,3),invf(lnp(100,2))),15)))
'''
n = 100 #nb de décimales à calculer
fc = fracont(racinek('2',n,12))
for k in range(20):
    rk = reduite(fc,k+1)
    print('Réduite d\'ordre '+str(k+1)+' de 2^(1/12): ',rk,'=',frac2dec(rk,20))
'''
