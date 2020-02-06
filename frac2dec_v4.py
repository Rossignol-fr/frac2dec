def explose_frac(r):
    """renvoie le numérateur et le dénominateur d'une fraction écrite sous
    forme de chaîne de caractères"""
    if '/' in r:
        pd = r.index('/') #position diviseur
        return int(r[:pd]), int(r[pd+1:])
    else:
        return int(r), 1

def explose_dec(r):
    """Transforme un décimal écrit sous forme de chaîne de caractères
    en liste sous la forme: [sgn,pE,pF,p, apx] où sgn est le signe (1 ou -1)
    pE la partie Entière, pF la partie décimale fixe, p la partie décimale
    périodique et apx un booléen indiquant si c'est une valeur approchée"""
    sgn = 1 # pour le signe
    apx = False
    if '~' in r: #Le nombre décimal est une valeur approchée
        r = r[1:] #on retire le flag
        apx = True
    if r[0] == '-':
        sgn = -1
        r = r[1:]
    if '.' not in r: # c'est un entier
        return sgn, r, '', '', apx
    else:
        pv = r.index('.') # donne la position du séparateur décimal
        pE = r[:pv]
    if '[' not in r:
        return sgn, pE, r[pv+1:], '', apx
    else:
        pf = r.index('[') # pour la partie décimale périodique
        return sgn, pE, r[pv+1:pf], r[pf+1:-1], apx

def frac2dec(r, n=-1):
    """Donne l'écriture décimale de la fraction p/q écrite
    sous la forme d'une chaîne de caractères. Si n est renseigné donne
    n chiffres après la virgule (par arrondi)"""
    p, q = explose_frac(r)
    sg, sgn = '', 1
    if p < 0:
        sg, sgn = '-', -1
        p = abs(p)
    r = p%q
    Lr = [r] #Liste des restes
    if r == 0: # C'est un entier !
        return frc((sgn*p//q,1))
    fs = sg + str(p//q) + '.'
    ls = len(fs)
    if n == 0: #on ne veut que la partie entière
        return fs[:-1]
    elif n == -1: #on veut tout le développement décimal
        j = n
        flg = False
        ic = 0
    else: # On ne veut que n chiffres après la virgule
        j = -n - 1 #On va chercher une décimale supplémentaire pour l'arrondi
        flg = True
        ic=1
    while j<0:
        j+=ic
        fs=fs+str(10*r//q)
        r=10*r%q
        if r in Lr:
            idx=Lr.index(r)
            fs=fs[:ls+idx]+'['+fs[ls+idx:]+']'
            break
        Lr.append(r)
        if r==0: # c'est un décimal !
            break
    if flg:
    #si on est là c'est une valeur approchée avec au moins une décimale
        fs = arrondi(fs,n)
        fs = chr(126) + fs #On place un tilde (~) devant le résultat
    return fs

#exemples: frac2dec('1/17') renvoie '0.[0588235294117647]'
# frac2dec('1/30') renvoie '0.0[3]'
# frac2dec('1/256') renvoie '0.00390625'
# frac2dec('17/12',n=5) renvoie '~1.41667'
# frac2dec('12/17',n=5) renvoie '~0.70588'

def simpfrac(a,b):
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
        return sgn*b, a
    else:
        return sgn*a, b

def pgcd(a,b):
    """Algorithme d'Euclide sur les entiers naturels a et b
    où a>b et renvoie le dernier reste non nul"""
    while a%b!=0:
        a, b= b, a%b
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
        f = simpfrac(int(dd)-int(d), sgn*(10**t-1)*10**s)
    else:
        nd = len(pF)
        f = simpfrac(int(pE)*10**nd+int(pF), sgn*10**nd)
    return frc(f)

def frc(f):
    """transforme un tuple d'entiers en fraction (string)"""
    if f[1] == 1:
        return str(f[0])
    else:
        return str(f[0]) + '/' + str(f[1])

def prodf(r1,r2,simp=True):
    """renvoie le produit de deux fractions écrites avec le symbole /"""
    p1, q1 = explose_frac(r1)
    p2, q2 = explose_frac(r2)
    if simp:
        return frc(simpfrac(p1*p2, q1*q2))
    else:
        return frc((p1*p2, q1*q2))

def addf(r1, r2, simp=True):
    """renvoie la somme de deux fractions écrites avec le symbole /"""
    p1, q1 = explose_frac(r1)
    p2, q2 = explose_frac(r2)
    if p1*q2 + q1*p2 == 0:
        return '0'
    else:
        if simp:
            return frc(simpfrac(p1*q2 + q1*p2, q1*q2))
        else:
            return frc((p1*q2 + q1*p2, q1*q2))

def extracine(a,n):
    """Renvoie la fraction approchant la racine carrée de a par l'algo
    de Babylone avec n décimales où a est un entier positif"""
    if a == 0 or a == 1:
        return a
    else:
        t = len(str(a))
        if t%2 == 0:
            u = 3*10**(t//2 - 1)
        else:
            u = 10**(t//2)
        s = str(u)
        e = addf(prodf(s, s), str(-a))
        j = 0
        while quotient_10n(e, n) and j < 10:
            p, q = explose_frac(s)
            s = prodf(addf(s, prodf(str(a), str(q)+'/'+str(p))), '1/2')
            e = addf(prodf(s, s), str(-a))
            j += 1
        return s, e

def quotient_10n(s, n):
    """Détermine si |s| > 10^-n où s est une fraction (string)"""
    p, q = explose_frac(s)
    if abs(q//p) < 10**n:
        return True
    else:
        return False

def arrondi(sd, n):
    """Renvoie l'arrondi à n décimales du nombre décimal (string)"""
    sgn, pE, pF, p, apx = explose_dec(sd)
    while len(pF) < n+1 and len(p) > 0:
        pF+=p
    if len(pF) < n + 1:
        return sd
    else:
        sn = pE + pF[:n + 1]
        sg=''
        pv=len(pE)
        if sgn==-1:
            sg='-'
        if int(sn[-1])>=5:
            nc=len(sn)
            bz='' #Pour le bourrage de zéros qui seront perdus
            nc1=len(str(int(sn)))
            nc2=len(str(int(sn)+5))
            i=0
            if nc2>nc1:
                i=1
            if nc1<nc:
                bz='0'*(nc-nc1+i)
            sn=bz+str(int(sn)+5)
            return sg + sn[:pv+i] + '.' + sn[pv+i:-1]
        else:
            return sg + pE + '.' + pF[:n]

def approxE(n):
    p,fac,f=1,1,'1'
    while fac<10**n:
        f=addf(f,frc((1,fac)))
        p+=1
        fac = fac*p
    return f

def approxPi2s6(n,simp=True):
    """serie des inverse des carrés (converge très lentement)"""
    p,f=1,'1'
    while p*p<10**n:
        p+=1
        f=addf(f,frc((1,p*p)),simp)
    return f