#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#USEFUL PACKAGES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import os,sys,commands,tempfile,time
from hashlib import md5
import numpy as np

import commands

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ALIASES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mod=np.linalg.norm
exit=sys.exit
argv=sys.argv
system=os.system
sleep=time.sleep
tmpf=tempfile.NamedTemporaryFile
body=dict
plot=dict

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#GLOBAL VARIABLES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROOTDIR=os.path.dirname(os.path.abspath(argv[0]))
PACKNAME="MercuPy"
VERBOSE=False
SIMULATE=False
OUTDIR="output"

################################################################################
#PHYSICAL QUANTITIES
################################################################################
#DOES NOT CHANGE MANTISA (VALUES USED BY MERCURY6). ORBITAL UNITS
PI=np.pi #Pi

GCONST=6.67428E-20 #km/s^2 kg #USNO
MSUN=1.9884E30 #kg #USNO
MEARTH=MSUN/332946.0487 #kg #USNO
MMARS=MSUN/3.09870359E6 #kg #USNO
MJUP=MSUN/1.047348644E3 #kg #USNO
MSAT=MSUN/3.4979018E3 #kg #USNO
MURANUS=MSUN/2.290298E4 #kg #USNO
MNEPTUNE=MSUN/1.941226E4 #kg #USNO

RSUN=6.96e5 #km
REARTH=6.371E3 #km
MMOON=7.3477E22 #kg
RMOON=1737.1 #km
RJUP=71492 #km
AU=1.49597870700e8 #km #USNO
HOURS=3600 #secs
DAYS=24*HOURS #secs
MONTHS=30*DAYS #secs
YEARS=365.25*DAYS #secs
GMC3=1E3 #kg/m^3
AUPD=AU/DAYS #km/s

################################################################################
#CUSTOM UNITS
################################################################################
G=1
UL=1.0*AU #km
UM=1.0*MSUN #kg
UT=np.sqrt(UL*UL*UL/(GCONST*UM)) #secs
UV=UL/UT #km/s
URHO=UM/(UL*UL*UL*1E9) #kg/m^3

################################################################################
#MERCURY UNITS
################################################################################
ULMERC=1.0*AU #km
UMMERC=1.0*MSUN #kg
UTMERC=1.0*DAYS #secs
UVMERC=ULMERC/UTMERC #km/s
GMERC=GCONST/(ULMERC**3/(UMMERC*UTMERC**2))
URHOMERC=1E3 #kg/m^3

################################################################################
#DATA INDEXING
################################################################################
ELSINDEX={'x':1,'y':2,'z':3,'vx':4,'vy':5,'vz':6,'r':7,
          'a':8,'e':9,'i':10,'g':11,'n':12,'M':13,
          'q':14,'Q':15,'l':16,'f':17,
          's':18,'o':19,'m':20,'d':21,
          'R':22,'P':23,'E':24,'L':25
          }

################################################################################
#OTHER CONSTANTS
################################################################################
EJECT_DISTANCE=100 #AU

################################################################################
#ENUMERATORS
################################################################################
#==================================================
#CLASS OF BODY
#==================================================
CUSTOMOBJECT=0
ROCKYPLSIMAL,ICEPLSIMAL=range(1,2+1)
NORMALSTAR,WDSTAR,NEUTRONSTAR=range(10,12+1)
GASGIANT,ICEGIANT=range(20,21+1)
SOLIDROCKY,SOLIDIRON,SOLIDICE=range(30,32+1)
#==================================================
#TYPE OF BODIES
#==================================================
CENTRAL,BIG,SMALL=range(0,3)
COORDINATES=dict()
#==================================================
#UNITS
#==================================================
"""
ORBITAL: km,kg,s,km/s
CUSTOM: UL,UM,UT,UV
MERCURY: AU,MSUN,DAYS,AU/DAY
BODY: RC,MC,s,km/s
GRAVITATIONAL: RH,MC,s,km/s
"""
ORBITAL,CUSTOM,MERCURY,BODY,GRAVITATIONAL=range(0,5)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#SYSTEM ROUTINES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def error(string,code=1):
    """
    Manage errors
    """
    print string
    print "Exiting %s."%PACKNAME
    exit(code)

class dictobj(object):
    """
    Class that allows the conversion from dictionary to class-like
    objects
    """
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        self.__dict__.update(other.__dict__)
        return self

def loadconf(filename):
    """Load configuration file

    Parameters:
    ----------
    filename: string
       Filename with configuration values.

    Returns:
    -------
    conf: dictobj
       Object with attributes as variables in configuration file

    Examples:
    --------
    >> loadconf('input.conf')

    """
    d=dict()
    conf=dictobj()
    if os.path.lexists(filename):
        execfile(filename,{},d)
        conf+=dictobj(d)
        qfile=True
    else:
        error("Configuration file '%s' does not found."%filename)
    return conf

def System(cmd,out=False,sim=SIMULATE):
    """
    Execute a command
    """
    if VERBOSE or sim:print "CMD:\n\t%s"%cmd
    if sim:return ""
    if not out:
        system(cmd)
        output=""
    else:
        output=commands.getoutput(cmd)
    return output

#Change template file to newfile
def Template2File(tempfile,key,val,newfile):
    """
    print "tempfile:",tempfile
    print "newfile:",newfile
    print "Key:",key
    print "Value:",val
    """
    system("sed -e 's/\[\[%s\]\]/%s/gi' %s > %s"%(key,val,tempfile,newfile))

#Change all field in a file for the dictionary entries
def Temp2File(tempfile,dic,newfile):
    tmp1=tmpf()
    tmp2=tmpf()
    System("cp -rf %s %s"%(tempfile,tmp1.name))
    for key in dic.__dict__:
        val=dic.__dict__[key]
        Template2File(tmp1.name,key,val,tmp2.name)
        System("cp -rf %s %s"%(tmp2.name,tmp1.name))
    System("cp -rf %s %s"%(tmp2.name,newfile))

#md5sum
def md5sum(string):
    md5s=md5(string)
    return md5s.hexdigest()

#Range
def Range(val,end=None,num=1,cyclic=False):
    if end is None:
        l=np.array([val])
    else:
        l=np.linspace(val,end,num)
    if cyclic:
        return l[:-1]
    else:
        return l

#Range
def List(ls):
    return np.array(ls)

#Get absolute directory
def absdir(directory):
    pwd=os.getcwd()
    os.chdir(directory)
    absd=os.getcwd()
    os.chdir(pwd)
    return absd

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PHYSICAL ROUTINES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def MassRadiusRelationship(Mp,ObjectClass):
    ScalingLaw="PowerLaw"

    #PLANETESIMALS
    if ObjectClass==ROCKYPLSIMAL:
        #REFERENCE OBJECT: CERES
        M1=9.43E20
        R1=487
        kappa=0.25
        pass
    elif ObjectClass==ICEPLSIMAL:
        #REFERENCE OBJECT: HALLEY COMET
        M1=2.2E14
        R1=15
        kappa=0.25
    #STARS
    elif ObjectClass==NORMALSTAR:
        #REFERENCE OBJECT: SUN
        M1=MSUN
        R1=RSUN
        kappa=0.8
    elif ObjectClass==WDSTAR:
        #REFERENCE OBJECT: SIRIUS B
        M1=0.978*MSUN
        R1=0.0084*RSUN
        kappa=-1./3
    elif ObjectClass==NEUTRONSTAR:
        #REFERENCE OBJECT: PSR B1257+12
        M1=1.5*MSUN
        R1=2E-5*RSUN
        kappa=-1./3
    #GIANT PLANETS
    elif ObjectClass==GASGIANT:
        #REFERENCE OBJECT: JUPITER
        M1=MJUP
        R1=RJUP
        kappa=0.25
    elif ObjectClass==ICEGIANT:
        #REFERENCE OBJECT: NEPTUNE
        M1=1.0243E26 #kg
        R1=24764 #km
        kappa=0.25
    #SOLID PLANETS
    elif ObjectClass==SOLIDROCKY:
        ScalingLaw="ModifiedLaw"
        M1=6.41*MEARTH
        R1=3.19*REARTH
        k1=-0.209594
        k2=+0.079900
        kappa=+0.413000
    elif ObjectClass==SOLIDICE:
        ScalingLaw="ModifiedLaw"
        M1=7.63*MEARTH
        R1=4.42*REARTH
        k1=-0.209396
        k2=+0.080700
        kappa=+0.375000
    elif ObjectClass==SOLIDIRON:
        ScalingLaw="ModifiedLaw"
        M1=4.34*MEARTH
        R1=2.23*REARTH
        k1=-0.209490
        k2=+0.080400
        kappa=+0.394000

    #SIMPLE POWER-LAW
    if ScalingLaw=="PowerLaw":
        Ms=Mp/M1
        Rs=Ms**kappa
        Rp=Rs*R1
    #MODIFIED POWER-LAW (Seager et al. 2008)
    elif ScalingLaw=="ModifiedLaw":
        Ms=Mp/M1
        Rs=10**(k1+1./3*np.log10(Ms)-k2*Ms**kappa)
        Rp=Rs*R1

    rho=Mp/(4*np.pi*(1E3*Rp)**3/3) #kg/m3
    return Rp,rho

#OBJVERB=True
OBJVERB=False
def ObjectProperties(body,fbody,cbody,COORDINATES):
    """
    body: dictionary
       ObjectClass: See mercurypprc.py
       Mass,Radius,Density

    Modify properties according to ObjectClass and values of Radius
    and Density (when ObjectClass is CustomObject).

    All values are in (SI): kg,km,g/cc

    Using mass determine the rest of properties
    """
    if OBJVERB:print "*"*90
    if OBJVERB:print "Info for body %s:"%body["Code"],body["Name"],body["Type"],body["ObjectClass"]
    
    #==============================
    #BASIC PROPERTIES
    #==============================
    #CUSTOM OBJECTS
    if body["ObjectClass"]==CUSTOMOBJECT:
        if body["Density"]>0:
            body["Radius"]=(body["Mass"]/(4*np.pi*body["Density"]/3))**(1./3)
            return 0
        elif body["Radius"]>0:
            body["Density"]=body["Mass"]/(4*np.pi*body["Radius"]**3)
            return 0
        else:
            error("Density and radius of Custom Object are both null.")

    #PREDEFINED OBJECTS
    body["Radius"],body["Density"]=MassRadiusRelationship(body["Mass"],
                                                          body["ObjectClass"])

    #CONVERSION TO MERCURY UNITS
    body["Mass"]/=UMMERC
    body["Radius"]/=ULMERC
    body["Density"]/=URHOMERC

    #==================================================
    #PROPERTIES RELATIVE TO CENTRAL BODY
    #==================================================
    body["Mu"]=GMERC*fbody["Mass"]

    if OBJVERB:print "Basic properties (M,R,rho,Mu):",body["Mass"],body["Radius"],body["Density"],body["Mu"]

    #==============================
    #CONVERSION FACTORS
    #==============================
    units=body["Units"]
    btype=body["Type"]
    if units==MERCURY:
        CL=1.0;CT=1.0;CV=1.0
    elif units==ORBITAL:
        CL=1.0/ULMERC
        CT=1.0/UTMERC
        CV=1.0/UVMERC
    elif units==CUSTOM:
        CL=UL/ULMERC
        CT=UT/UTMERC
        CV=UV/UVMERC
    elif units==BODY:
        CL=body["Radius"]
        CT=1.0/UTMERC
        CV=1.0/UVMERC
    elif units==GRAVITATIONAL:
        CL=fbody["RH"]
        CT=1.0/UTMERC
        CV=1.0/UVMERC

    if OBJVERB:print "System of coordinates:",COORDINATES[btype]
    if OBJVERB:print "Units conversion (CL,CT,CV):",CL,CT,CV

    #==============================
    #STATE VECTOR -> CSPICE FORMAT
    #==============================
    stvec=np.array([float(st) for st in body["State"].split()])
    qtransform=False
    if "Ast" in COORDINATES[btype]:
        stvec[0]*=CL
        q=stvec[0]*(1-stvec[1])
        state="%e "%q
        state+="%e %e %e %e %e "%tuple(stvec[1:6])
        state+="%e %e %e"%(body["Ep"],body["Mu"],body["Ep"])
        qtransform=True

    elif "Com" in COORDINATES[btype]:
        stvec[0]*=CL
        stvec[-1]*=CT
        state="%e "%stvec[0]
        state+="%e %e %e %e %e "%tuple(stvec[1:6])
        state+="%e %e %e"%(body["Ep"],body["Mu"],body["Ep"])
        qtransform=True

    elif "Car" in COORDINATES[btype]:
        state=body["State"]
        #print stvec
        stvec[0:3]*=CL
        stvec[3:6]*=CV
        statec="%e %e %e %e %e %e"%tuple(stvec)
        qtransform=False

    if OBJVERB:print "State (Mercury units):",stvec
    if 'Cartesian' in body:
        if OBJVERB:print "Cartesian:",body["Cartesian"]
        pass
    
    if OBJVERB:print "CSPICE transformation:",qtransform
    #==================================================
    #CONVERT ALL POSITIONS TO CENTRAL BODY POSITIONS
    #==================================================
    if OBJVERB:print "State for CSPICE program:",state
    body["Units"]=MERCURY
    body["Frame"]=cbody["Code"]

    #CSPICE CONVERSION FROM ELEMENTS TO STATE
    if OBJVERB:print "State previous to CSPICE:",state
    if qtransform:
        statec=System("bin/elem2state %s"%state,out=True,sim=False)

    if OBJVERB:print "State coordinates:",statec
    if OBJVERB:print "State from CSPICE:",statec
   
    #MOVE COORDINATES TO CENTRAL BODY
    body["Cartesian"]=np.array([float(st) for st in statec.split()])
    if OBJVERB:print "Cartesian respect fbody:",body["Cartesian"]
    if OBJVERB:print "Cartesian of fbody:",fbody["Cartesian"]
    body["Cartesian"]+=fbody["Cartesian"]
    if OBJVERB:print "Cartesian respect central body:",body["Cartesian"]

    #STORE FINAL STATE RESPECT TO CENTRAL BODY
    state="%e %e %e %e %e %e"%tuple(body["Cartesian"])
    body["State"]=state
    if OBJVERB:print "Final State:",body["State"]

    #==================================================
    #COMPUTE THE HILL RADIUS RESPECT TO FRAME BODY
    #==================================================
    relpos=body["Cartesian"]-fbody["Cartesian"]
    d=np.sqrt((relpos[0:2]**2).sum())
    if OBJVERB:print "Distance to frame body:",d
    if body["Type"]!=CENTRAL:
        body["RH"]=d*(body["Mass"]/(3*fbody["Mass"]))**(1./3)
    else:
        body["RH"]=EJECT_DISTANCE
    if OBJVERB:
        print "Hill distance respect frame body:",body["RH"]
        print "*"*90
    return 0

def ShowBody(body,fhl=sys.stdout):
    """
    Show body properties
    """
    print>>fhl, "Object %d (Type,Obj.Class): Code %s, Name %s (%d,%d)"%(body["Id"]+1,body["Code"],body["Name"],body["Type"],body["ObjectClass"])
    print>>fhl, "\tMass (UM = %e kg): %e UM, %e kg"%(UMMERC,body["Mass"],body["Mass"]*UMMERC)
    print>>fhl, "\tRadius (UL = %e km): %e UL, %e km"%(ULMERC,body["Radius"],body["Radius"]*ULMERC)
    print>>fhl, "\tDensity (URHO = %e g/cc): %e URHO, %e kg/m3"%(URHOMERC,body["Density"],body["Density"]*URHOMERC)
    print>>fhl, "\tClose Encounter (UL = %e km): %e RH, %e UL, %e km"%(ULMERC,body["CloseEncounter"],body["CloseEncounter"]*body["RH"],body["CloseEncounter"]*body["RH"]*ULMERC)
    print>>fhl, "\tFrame object : %s"%body["Frame"]
    print>>fhl, "\tMu : %e"%body["Mu"]
    print>>fhl, "\tHill Radius (UL = %e km): %e UL, %e km"%(ULMERC,body["RH"],body["RH"]*ULMERC)
    print>>fhl, "\tOriginal coordinates : %s"%COORDINATES[body["Type"]]
    state=np.array([st for st in body["Cartesian"]])
    states="%+e UL %+e UL %+e UL %+e UV %+e UV %+e UV"%tuple(state)
    print>>fhl, "\tCartesian position (UL=%e km,UT=%e s,UV=%e km/s):\n\t\t"%(ULMERC,UTMERC,UVMERC),states
    state[0:3]*=ULMERC
    state[3:6]*=UVMERC
    states="%+e km %+e km %+e km %+e km/s %+e km/s %+e km/s"%tuple(state)
    print>>fhl, "\t\t",states
    print>>fhl

def rotate3D(R,phi,theta):
    from scipy import dot
    X,Y,Z=0,1,2
    #ROTATION MATRIX
    p=phi*PI/180
    t=theta*PI/180
    c=np.cos
    s=np.sin
    R3=np.array(
        [
            [c(p),s(p),0],
            [-s(p)*c(t),c(p)*c(t),s(t)],
            [s(p)*s(t),-c(p)*s(t),c(t)]
         ]
        );
    Xv=[]
    Yv=[]
    Zv=[]
    for r in R:
        rv=dot(R3,r)
        Xv+=[rv[X]]
        Yv+=[rv[Y]]
        Zv+=[rv[Z]]

    return Xv,Yv,Zv

