#!/usr/bin/env python

#
# TODO:
#   o Support array(swizzle) much more cleaner manner.
#   o and more...
#

import os
import os.path
import subprocess
import sys
import tempfile
from string import Template

from sexps import *

gSlotToN = {
    'x' : 0
  , 'y' : 1
  , 'z' : 2
  , 'w' : 3 
}

gNToSlot = ['x', 'y', 'z', 'w']

gIndentLevel = 0
gVarCount = 0

# Symbol table(as stack) for tempolrary variables
# depth 1 = global scope
# depth 2 = function scope
gSymbolTable        = []

# @todo { Create context to save these variables. }
gUniformInputs      = {}
gUniformCount       = 0

gVaryingInputs      = {}
gVaryingInputCount  = 0

gVaryingOutputs     = []

gStructDefs         = {}

def DefineStruct(struct):
    gStructDefs[struct['name']] = struct

def FindStruct(name):
    if isinstance(name, list) and name[0] == 'array':
        # maybe array of struct
        # ['array', 'Sphere__0x7fdac1c11eb0', '3']
        p = name[1].split('__')
        if len(p) > 1:
            sname = p[0] + '__' 
        else:
            sname = p[0] 
        if gStructDefs.has_key(sname):
            return gStructDefs[sname]
        raise
    elif gStructDefs.has_key(name):
        return gStructDefs[name]

    return False

def GetGLTypeSize(basety, n):
    # 4 = sizeof(float)
    # 8 = sizeof(double)
    if basety == "vec":
        return 4 * n
    if basety == "dvec":
        return 8 * n
    elif basety == "ivec":
        return 4 * n
    elif basety == "bvec":
        return 1 * n    # byte
    elif basety == "mat":
        return 4 * n * n
    
    # Unknown type
    assert 0
    
def AddSymbol(varname, ty, n, quals):
    assert len(gSymbolTable) > 0

    # Add symbol to last stack of symbol table
    gSymbolTable[-1][varname] = {'type' : (ty, n), 'quals' : quals, 'name' : varname}

def GetSymbol(varname):
    assert len(gSymbolTable) > 0

    # search internal scope first
    n = len(gSymbolTable) - 1 
    while n >= 0:

        if varname in gSymbolTable[n]:
            return gSymbolTable[n][varname]

        n = n - 1


    return False

def GetScopeLevel():
    return len(gSymbolTable)
    

def IsArray(s):
    return isinstance(s, (list, tuple))

def IsSamplerTy(tyname):
    samplers = [
        "sampler2D"
      , "sampler3D"
    ]

    if tyname in samplers:
        return True

    return False

def ParseTy(tyname):
    """
    Parse type into (elementType, num). sampler has special treatment
      example. vec4      -> ('vec'  , 4)
               float     -> ('float', 1)
               sampler2D -> ('sampler2D', 1)
    """

    if IsSamplerTy(tyname):
        return (tyname, 1)

    # array type is TODO
    if IsArray(tyname):
        return (tyname, 1)

    p = re.compile("([a-zA-z]+)([0-9]*)")

    groups = p.match(tyname).groups()
    if len(groups) == 1:
        return (groups[0], 1)
    else:
        if groups[1] == '':
            return (groups[0], 1)
        else:
            return (groups[0], int(groups[1]))


def IsBuiltinTextureFunction(fname):
    builtins = [
        "texture2D"
      , "texture3D"
    ]

    if fname in builtins:
        return True

    return False

def IsBuiltinTraceFunction(fname):
    builtins = [
        "trace"                 # LSGL ext
      , "raydepth"              # LSGL ext
      , "rayattrib"             # LSGL ext
      , "rayoption"             # LSGL ext
      , "isectinfo"             # LSGL ext
      , "camerainfo"            # LSGL ext
      , "numIntersects"         # LSGL ext
      , "queryIntersect"        # LSGL ext
    ]

    if fname in builtins:
        return True

    return False

def IsBuiltinRandomFunction(fname):
    builtins = [
        "random"               # LSGL ext
    ]

    if fname in builtins:
        return True

    return False

def IsBuiltinFunction(fname):
    builtins = [
        "radians"
      , "degrees"
      , "sin"
      , "cos"
      , "tan"
      , "asin"
      , "acos"
      , "atan"
      , "atan2"
      , "pow"
      , "exp"
      , "log"
      , "exp2"
      , "log2"
      , "sqrt"
      , "inversesqrt"
      , "abs"
      , "sign"
      , "floor"
      , "ceil"
      , "fract"
      , "mod"
      , "min"
      , "max"
      , "clamp"
      , "mix"
      , "step"
      , "smoothstep"
      , "length"
      , "distance"
      , "dot"
      , "cross"
      , "normalize"
      , "faceforward"
      , "reflect"
      , "refract"
      , "matrixCompMult"
      , "lessThan"
      , "lessThanEqual"
      , "greaterThan"
      , "greaterThanEqual"
      , "equal"
      , "notEqual"
      , "any"
      , "all"
      , "not"
      , "texture2D"
      , "texture3D"
      , "trace"                 # LSGL ext
      , "raydepth"              # LSGL ext
      , "rayattrib"             # LSGL ext
      , "rayoption"             # LSGL ext
      , "isectinfo"             # LSGL ext
      , "camerainfo"            # LSGL ext
      , "random"                # LSGL ext
      , "numIntersects"         # LSGL ext
      , "queryIntersect"        # LSGL ext
    ]

    if fname in builtins:
        return True

    return False

def GetBuiltinType(varname):
    builtins = {
        "gl_FragCoord" : ("vec", 4)
      , "gl_FrontFacing" : ("bool", 1)
      , "gl_FragColor" : ("vec", 4)
      , "gl_PointCoord" : ("vec", 2)
      , "gl_MaxVertexAttribs" : ("int", 1)
      , "gl_MaxVertexUniformVectors" : ("int", 1)
      , "gl_MaxVaryingVectors" : ("int", 1)
      , "gl_MaxVertexTextureImageUnits" : ("int", 1)
      , "gl_MaxCombinedTextureImageUnits" : ("int", 1)
      , "gl_MaxTextureImageUnits" : ("int", 1)
      , "gl_MaxFragmentUniformVectors" : ("int", 1)
      , "gl_DepthRangeParameters" : ("todo", 1)
      , "gl_DepthRange" : ("todo", -1)
      , "gl_MaxDrawBuffers" : ("int", 1)
      , "gl_FragData" : ("todo", -1)
    }

    if varname in builtins:
        return builtins[varname]

    return False

def IsBuiltinVariable(varname):
    builtins = [
        "gl_FragCoord"
      , "gl_FrontFacing"
      , "gl_FragColor"
      , "gl_PointCoord"
      , "gl_MaxVertexAttribs"
      , "gl_MaxVertexUniformVectors"
      , "gl_MaxVaryingVectors"
      , "gl_MaxVertexTextureImageUnits"
      , "gl_MaxCombinedTextureImageUnits"
      , "gl_MaxTextureImageUnits"
      , "gl_MaxFragmentUniformVectors"
      , "gl_DepthRangeParameters"
      , "gl_DepthRange"
      , "gl_MaxDrawBuffers"
      , "gl_FragData"
    ]

    if varname in builtins:
        return True

    return False

def AddUniformInput(varname, ty, n, quals):

    # Assign unique index(from 0)
    global gUniformCount
    gUniformInputs[varname] = {'type' : (ty, n), 'quals' : quals, 'name' : varname, 'index' : gUniformCount}
    gUniformCount += 1

def IsUniform(varname):
    global gUniformInputs
    if varname in gUniformInputs:
        return True

    return False

def GetUniform(varname):
    global gUniformInputs

    return gUniformInputs[varname]

def IsVaryingInput(varname):
    global gVaryingInputs
    if varname in gVaryingInputs:
        return True

    return False

def GetVaryingInput(varname):
    global gVaryingInputs

    return gVaryingInputs[varname]

def IsTexture(varname):

    if IsUniform(varname):
        uniform = GetUniform(varname)

        if IsSamplerTy(uniform['type'][0]):
            return True

    return False

def IsTemporary(varname):
    if IsTexture(varname) or IsVaryingInput(varname) or IsUniform(varname):
        return False

    # Might be temporary variable
    return True

def AddVaryingInput(varname, ty, n, quals):
    # Assign unique index(from 0)
    global gVaryingInputCount
    gVaryingInputs[varname] = {'type' : (ty, n), 'quals' : quals, 'name' : varname, 'index' : gVaryingInputCount}
    gVaryingInputCount += 1

def GetTypeCastString(varname):

    if IsVaryingInput(varname):
        varying = GetVaryingInput(varname)

        if IsVectorType(varying['type'][0]):
            return "(%s%d *)" % (varying['type'][0], varying['type'][1])
        elif IsMatrixType(varying['type'][0]):
            return "(%s%d *)" % (varying['type'][0], varying['type'][1])
        else:
            return "(%s *)" % (varying['type'][0])

    elif IsUniform(varname):
        uniform = GetUniform(varname)

        if IsTexture(varname):
            # No typecast required
            return ""
        else:
            if IsVectorType(uniform['type'][0]):
                return "(%s%d *)" % (uniform['type'][0], uniform['type'][1])
            elif IsMatrixType(uniform['type'][0]):
                return "(%s%d *)" % (uniform['type'][0], uniform['type'][1])
            else:
                return "(%s *)" % (uniform['type'][0])

    # No typecast required
    return ""

def GetTypeOfSymbol(varname):

    if IsVaryingInput(varname):
        varying = GetVaryingInput(varname)

        return varying['type']

    elif IsUniform(varname):
        uniform = GetUniform(varname)

        return uniform['type']

    elif IsBuiltinVariable(varname):
        return GetBuiltinType(varname)

    elif IsTemporary(varname):
        temp = GetSymbol(varname)

        # print varname
        assert temp is not False
        return temp['type']

    assert 0    # Unknown symbol


def IncrementIndent():
    global gIndentLevel
    gIndentLevel += 1

def DecrementIndent():
    global gIndentLevel
    gIndentLevel -= 1
    if gIndentLevel < 0:
        gIndentLevel = 0

def Indent():
    global gIndentLevel
    s = ""
    for i in range(gIndentLevel):
        s += "  "

    return s
        
def NewTempVar():
    global gVarCount

    s = "tmpvar_%d" % gVarCount

    gVarCount += 1

    return s

def IsVectorType(s):
    tys = [ 'vec', 'ivec', 'bvec', 'dvec' ]
    if s in tys:
        return True
    
    return False

def IsMatrixType(s):
    tys = [ 'mat' ]
    if s in tys:
        return True
    
    return False

def baseType(ty):

    if ty == "vec":
        return "float"
    elif ty == "ivec":
        return "int";
    elif ty == "bvec":
        return "bool"
    elif ty == "mat":
        return "float"
    else:
        return ty

def parseValue(lst):

    (ty, n) = ParseTy(lst[0])
    values = lst[1]

    return (ty, n, values)

def renameVariable(varname):
    #if varname == "gl_FragColor":
    #    return "(fragment->fragColor)"

    # Rewrite builtin variables
    if varname == "gl_FragColor":
        return "(__fragment->fragColor)"
    elif varname == "gl_FragCoord":
        return "(__fragment->fragCoord)"
    elif varname == "Normal":
        return "(__fragment->normal)"
    elif varname == "Eye":
        return "(@TODO:Eye)"

    # Rewrite uniform variables
    if IsUniform(varname):
        uniform = GetUniform(varname)

        if IsTexture(varname):
            return "(__state->textures[%d])" % uniform['index']
        else:
            return "(__state->uniforms[%d])" % uniform['index']

    # Rewrite varying variables
    if IsVaryingInput(varname):
        varying = GetVaryingInput(varname)

        return "(__state->varyings[%d])" % varying['index']

    # Might be local variable. No rewrite required.
    return varname

class VarDecl:
    def __init__(self, varname, ty, n):
        self.varname = varname
        self.ty      = ty
        self.n       = n

    def __str__(self):
        return self.varname

    def getDeclareString(self):
        s = ""

        if isinstance(self.ty, list) and self.ty[0] == 'array':
            # array def
            # ['array', 'Sphere__0x7fe31ac11eb0', '3']
            p = self.ty[1].split('__')
            if len(p) > 1:
                sname = p[0] + '__'
            else:
                sname = p[0]
            s += "%s %s[%s];\n" % (sname, self.varname, self.ty[2])
        elif IsVectorType(self.ty):
            s += "%s%d %s;\n" % (self.ty, self.n, self.varname)
        elif IsMatrixType(self.ty):
            s += "%s%d %s;\n" % (self.ty, self.n, self.varname)
        else:
            s += "%s %s;\n" % (self.ty, self.varname)

        return s

    def getIntermExprString(self):
        return ""

    def getExprString(self, slot, i):
        return self.varname + ("[%d]" % gSlotToN[slot])
            

class VarRef:
    def __init__(self, varname):
        self.varname = renameVariable(varname)
        self.orgname = varname

    def __str__(self):
        return self.varname

    def getDeclareString(self):
        return ""

    def getIntermExprString(self):
        return ""

    def getCExpr(self):

        prefix = ""
        isInput = False

        if IsUniform(self.orgname) or IsVaryingInput(self.orgname):
            # If this variable is uniform/varying variable, add a prefix.
            prefix = ".data"
            isInput = True

        (ty, n) = GetTypeOfSymbol(self.orgname)
        tycast = GetTypeCastString(self.orgname)

        if isInput:
            return "/*var_ref*/" + "(*" + tycast + "(" + self.varname + prefix + "))"

        else:

            postfix = "&"

            return "/*var_ref*/" + "(*" + tycast + postfix + "(" + self.varname + prefix + "))"

    def getExprString(self, slot, i):
        prefix = ""
        isInput = False

        if IsUniform(self.orgname) or IsVaryingInput(self.orgname):
            # If this variable is uniform variable, add a prefix.
            prefix = ".data"
            isInput = True

        else:
            # If this variable is temporary variable, add some prefix depending on its type.
            var = GetSymbol(self.orgname)
            if var is not False:
                if IsVectorType(var['type'][0]):
                    prefix = ".v"

        (ty, n) = GetTypeOfSymbol(self.orgname)
        tycast = GetTypeCastString(self.orgname)

        if isInput:
            # print (ty, n)
            if IsVectorType(ty):
                return "/*var_ref(vec)*/" + "(*" + tycast + "(" + self.varname + prefix + ")).v[%d]" % gSlotToN[slot]
            else:
                return "/*var_ref*/" + "(*" + tycast + "(" + self.varname + prefix + "))"

        else:

            postfix = "&"

            if IsVectorType(ty):
                return "/*var_ref*/" + tycast + self.varname + prefix + ("[%d]" % gSlotToN[slot])
            else:
                return "/*var_ref*/" + "(*" + tycast + postfix + "(" + self.varname + prefix + "))"
            
class RecordRef:
    def __init__(self, varname, membername):

        if (len(varname) == 2 and varname[0] == 'var_ref'):
            # should be ['var_ref', 'str']
            self.var = VarRef(varname[1])
            self.recordname = varname[1]
            self.membername = membername
            self.is_array = False
        elif (len(varname) == 3 and varname[0] == 'array_ref'):
            # ['array_ref', ['var_ref', 'sphere'], ['constant', 'int', ['0']]]
            self.var = ArrayRef(varname[1], varname[2])
            self.recordname = varname[1][1]
            self.membername = membername
            self.is_array = True
        else:
            raise
    def __str__(self):
        return str(self.var) + "." + self.membername

    def getDeclareString(self):
        return ""

    def getIntermExprString(self):
        return ""

    def getCExpr(self):

        # Look up struct definition
        sym = GetSymbol(self.recordname)
        assert sym is not False

        (sty, sn) = sym['type']

        struct = FindStruct(sty)
        assert struct is not False


        # Look up type of member variable
        member = struct['members'][self.membername]
        assert member is not False

        (ty, n) = ParseTy(member['ty'])

        return "/*record_ref*/" + "(" + str(self.var) + "." + self.membername + ")" 

    def getExprString(self, slot, i):

        # Look up struct definition
        sym = GetSymbol(self.recordname)
        assert sym is not False

        (sty, sn) = sym['type']

        # print((sty, sn))
        struct = FindStruct(sty)
        assert struct is not False


        # Look up type of member variable
        member = struct['members'][self.membername]
        assert member is not False

        (ty, n) = ParseTy(member['ty'])

        # return "/*record_ref*/%s" % struct['name']

        prefix = ""

        if IsVectorType(ty):
            prefix = ".v"

        if IsVectorType(ty):
            return "/*record_ref*/" + "(" + str(self.var) + "." + self.membername + ")" + prefix + ("[%d]" % gSlotToN[slot])
        else:
            return "/*record_ref*/" + "(" + str(self.var) + "." + self.membername + ")"


class ArrayRef:
    def __init__(self, varname, indexname):
        assert len(varname) == 2;   # should be ['var_ref', 'str']
        assert varname[0] == 'var_ref';

        assert len(indexname) == 3;   # should be ['constant', 'int', [N]]

        self.var = VarRef(varname[1])
        self.recordname = varname[1]
        self.indexname = indexname[2][0]

    def __str__(self):
        return self.recordname + "[" + self.indexname + "]"

    def getDeclareString(self):
        return ""

    def getIntermExprString(self):
        return ""

    def getCExpr(self):

        # Look up struct definition
        sym = GetSymbol(self.recordname)
        assert sym is not False

        #(sty, sn) = sym['type']

        #structinfo = FindStruct(sty)
        #assert structinfo is not False

        return "/*array_ref*/" + self.recordname + "[%s]" % self.indexname

    def getExprString(self, slot, i):

        if IsUniform(self.recordname):
            # Assume mat3, mat4 or elase

            prefix = ".data"

            (ty, n) = GetTypeOfSymbol(self.recordname)
            tycast = GetTypeCastString(self.recordname)
            name = renameVariable(self.recordname)

            if IsMatrixType(ty):
                return "/*array_ref*/" + "(*" + tycast + "(" + name + prefix + (")).v[%s][%d]" % (self.indexname, gSlotToN[slot]))
            else:
                print (ty, n)
                raise
                #return "/*array_ref*/" + "(" + self.recordname + "[" + self.indexname + "]" + ")"

        else:

            # Look up struct definition
            sym = GetSymbol(self.recordname)
            #print sym, self.recordname
            assert sym is not False

            (ty, n) = sym['type']

            prefix = ""

            if IsVectorType(ty):
                prefix = ".v"
            elif IsMatrixType(ty):
                prefix = ".v"

            if IsVectorType(ty):
                return "/*array_ref*/" + "(" + self.recordname + "[" + self.indexname + "]" + prefix + ("[%d]" % gSlotToN[slot]) + ")"
            if IsMatrixType(ty):
                return "/*array_ref*/" + "(" + self.recordname + prefix + ("[%s][%d]" % (self.indexname, gSlotToN[slot])) + ")"
            else:
                return "/*array_ref*/" + "(" + self.recordname + "[" + self.indexname + "]" + ")"

class Constant:
    def __init__(self, ty, n, values):
        self.ty = ty
        self.n  = n
        i = 0;

        self.values = []
        # Consider '10', 'e', '-10' case
        while i < len(values):
            if i + 2 < len(values) and (values[i+1] == 'e' or value[i+1] == 'E'):     
                val = values[i] + values[i+1] + values[i+2]
                self.values.append(val)
                i = i + 3
            else:
                self.values.append(values[i])
                i = i + 1

    def __str__(self):
        return "(%s, %d, %s)" % (self.ty, self.n, self.values)

    def getDeclareString(self):
        return ""

    def getIntermExprString(self):
        return ""

    def getCExpr(self):
        s = ""

        if (self.n == 1):
            s += "(" + str(self.values[0]) + ")"
        else:
            if IsVectorType(self.ty):
                # for const value, no swizzle required.
                # n = gSlotToN[slot]
                # s = str(self.values[i])
                s += "__make_vec%d(" % self.n
                s += ','.join(map(str, self.values))
                s += ")"
            elif IsMatrixType(self.ty):
                # for const value, no swizzle required.
                # n = gSlotToN[slot]
                # s = str(self.values[i])
                s += "__make_mat%d(" % self.n
                s += ','.join(map(str, self.values))
                s += ")"
            else:
                raise

        return s

    def getExprString(self, slot, i):
        s = ""

        if (self.n == 1):
            s = "(" + str(self.values[0]) + ")"
        else:
            # for const value, no swizzle required.
            # n = gSlotToN[slot]
            s = "(" + str(self.values[i]) + ")"

        return s


class Swizzle:
    def __init__(self, slots, values):
        self.slots = slots
        self.values = values;

    def __str__(self):
        #return "TODO(swizzle: %s, %s)" % (self.slots, str(self.values))
        return ""

    def getDeclareString(self):
        s = ""
        s += self.values.getDeclareString()

        return s

    def getIntermExprString(self):

        s = ""
        s += self.values.getIntermExprString()

        #for i in range(self.n):
        #    slot = gNToSlot[i]
        #    s += Indent()
        #    s += self.getExprString(slot, i)
        #    s += " = "
        #    s += self.lhs.getExprString(slot, i)
        #    s += " " + self.op + " "
        #    s += self.rhs.getExprString(slot, i)
        #    s += ";\n"
        return s

    def getCExpr(self):
        s = ""
        s += "__swizzle("
        for (n, i) in enumerate(self.slots):
            s += "%d" % gSlotToN[i]
            if n != len(self.slots):
                s += ", "
        s += self.values.getCExpr()
        s += ")"
        return s

    def getExprString(self, slot, i):
        s = ""
        nn = len(self.slots)
        if (nn == 1):
            s = "/*swiz*/"
            s += self.values.getExprString(self.slots, 0)
            #s += "[%d]" % gSlotToN[self.slots]
            assert isinstance(s, str)
        else:
            #n = gSlotToN[slot]
            ss = self.slots[i]
            s = self.values.getExprString(ss, i)
            assert isinstance(s, str)

        return s

class Break:
    def __init__(self):
        pass

    def __str__(self):
        return "break"

    def getDeclareString(self):
        return ""

    def getIntermExprString(self):
        return ""

    def getCExpr(self):
        return "break"

    def getExprString(self, slot, i):
        return "break"

class Assign:
    def __init__(self, slots, lhs, rhs):
        self.slots = slots
        self.lhs   = lhs;
        self.rhs   = rhs;

    def __str__(self):
        #return "TODO(swizzle: %s, %s)" % (self.slots, str(self.values))
        return ""

    def getDeclareString(self):
        return ""

    def getIntermExprString(self):
        return ""
        #s = ""
        #s += self.values.getIntermExprString()

        #for (n, slot) in enumerate(self.slots):
        #    #slot = gNToSlot[i]
        #    s += Indent() + self.values.getExprString(slot, n) + ";\n"
        #return s

    def getCExpr(self):
        # todo
        return ""

    def getExprString(self, slot, i):
        s = ""
        nn = len(self.slots)
        if (nn == 1):
            s = str(self.values)
            s += "[%d]" % i
            assert isinstance(s, str)
        else:
            #n = gSlotToN[slot]
            ss = self.slots[i]
            s = self.values.getExprString(ss, i)
            assert isinstance(s, str)

        return s

class BinaryExpression:
    def __init__(self, ty, n, op, lhs, rhs):
        self.ty     = ty
        self.n      = n
        self.op     = op
        self.lhs    = lhs
        self.rhs    = rhs
        self.dst    = NewTempVar()

    def __str__(self):
        return "TODO(expr)"

    def getDeclareString(self):

        s = ""
        s += self.lhs.getDeclareString()
        s += Indent() + self.rhs.getDeclareString()

        baseTy = baseType(self.ty)
        if self.n == 1:
            s += "%s %s;\n" % (baseTy, self.dst)
        else:
            s += "%s %s[%d];\n" % (baseTy, self.dst, self.n)

        return s

    def getIntermExprString(self):
        s = ""
        s += self.lhs.getIntermExprString()
        s += self.rhs.getIntermExprString()

        if self.ty == 'vec' and self.n >= 2: # mat and vec?
            ops = {
                '/' : '__div'
              , '+' : '__add'
              , '-' : '__sub'
              , '*' : '__mul'
              , '>' : '__gt'
              , '>=' : '__ge'
              , '<' : '__lt'
              , '<=' : '__le'
              , '==' : '__eq'
              , '&&' : '__and'
              , '||' : '__or'
              , '!' : '__not'
              , 'any_nequal' : '__any_neq'
              , 'all_equal' : '__all_eq'
            }

            func = ops[self.op]
            assert func is not False

            s = ""
            s += func + "("
            s += self.lhs.getCExpr()
            s += ", "
            s += self.rhs.getCExpr()
            s += ");\n"

        else:
            for i in range(self.n):
                slot = gNToSlot[i]
                s += Indent()
                s += self.getExprString(slot, i)
                s += " = "
                s += self.lhs.getExprString(slot, i)
                s += " " + self.op + " "
                s += self.rhs.getExprString(slot, i)
                s += ";\n"

        return s

    def getCExpr(self):

        ops = {
            '/' : '__div'
          , '+' : '__add'
          , '-' : '__sub'
          , '*' : '__mul'
          , '>' : '__gt'
          , '>=' : '__ge'
          , '<' : '__lt'
          , '<=' : '__le'
          , '==' : '__eq'
          , '&&' : '__and'
          , '||' : '__or'
          , '!' : '__not'
          , 'any_nequal' : '__any_neq'
          , 'all_equal' : '__all_eq'
        }

        func = ops[self.op]
        assert func is not False

        s = ""
        s += func + "("
        s += self.lhs.getCExpr()
        s += ", "
        s += self.rhs.getCExpr()
        s += ")"

        return s

    def getExprString(self, slot, i):
        s = ""
        if self.n == 1:
            s = self.dst
        else:
            s = "%s[%d]" % (self.dst, gSlotToN[slot])
        assert isinstance(s, str)
        return s

    def dst(self):
        return self.dst

class UnaryExpression:
    def __init__(self, ty, n, op, src):
        self.ty     = ty
        self.n      = n
        self.op     = op
        self.src    = src
        self.dst    = NewTempVar()

    def __str__(self):
        return "TODO(uexpr)"

    def getDeclareString(self):
        s = ""
        s += self.src.getDeclareString()

        baseTy = baseType(self.ty)
        if self.n == 1:
            s += "%s %s;\n" % (baseTy, self.dst)
        else:
            s += "%s %s[%d];\n" % (baseTy, self.dst, self.n)

        return s

    def getIntermExprString(self):
        # print "UnaryInterm\n"

        ops = {
            'neg' : '-('
          , 'rcp' : '__rcp('
          , 'i2f' : '(float)('
          , 'f2i' : '(int)('
          , 'b2f' : '(float)('
          , '!' : '!('
        }

        opExpr = ops[self.op]
        assert opExpr is not False

        s = ""
        s += self.src.getIntermExprString()

        for i in range(self.n):
            slot = gNToSlot[i]
            s += Indent()
            s += self.getExprString(slot, i)
            s += " = "
            s += ("(%s(" % opExpr) + self.src.getExprString(slot, i) + ")))"
            s += ";\n"
        return s

    def getCExpr(self):

        ops = {
            'neg' : '__neg'
          , 'rcp' : '__rcp'
          , 'rsq' : '__rsq'
          , 'i2f' : '__i2f'
          , 'f2i' : '__f2i'
          , 'b2f' : '__b2f'
          , '!' : '__not'
        }

        func = ops[self.op]
        assert func is not False

        s = ""
        s += func + "("
        s += self.src.getCExpr()
        s += ")"

        return s

    def getExprString(self, slot, i):

        s = ""
        if self.n == 1:
            s = "%s" % (self.dst)
        else:
            s = "%s[%d]" % (self.dst, gSlotToN[slot])
        assert isinstance(s, str)
        return s

    def dst(self):
        return self.dst

def constructExpr(expr):

    name = expr[0]

    if name == "var_ref":
        return VarRef(expr[1])
    elif name == "record_ref":
        return RecordRef(expr[1], expr[2])
    elif name == "constant":
        (ty, n, values) = parseValue(expr[1:])
        return Constant(ty, n, values)
    elif name == "assign":
        slots = expr[1]
        lhs   = constructExpr(expr[2])
        rhs   = constructExpr(expr[3])
        return Assign(slots, lhs, rhs)
    elif name == "swiz":
        slots = expr[1]
        e = constructExpr(expr[2])
        return Swizzle(slots, e)
    elif name == "expression":
        (ty, n) = ParseTy(expr[1])

        op = expr[2]
        lhs = constructExpr(expr[3])
        if len(expr) > 4:
            rhs = constructExpr(expr[4])
            return BinaryExpression(ty, n, op, lhs, rhs)
        else:
            # Unary expression
            return UnaryExpression(ty, n, op, lhs)

    elif name == "declare":
        quals = expr[1]
        ty = expr[2]
        offt = 3;

        (ty, n) = ParseTy(ty)

        varname = expr[offt]

        return VarDecl(varname, ty, n)

    elif name == "call":
        # print expr
        assert 0

    elif expr == "break":   # not a name
        return Break()
    elif name == "array_ref":
        return ArrayRef(expr[1], expr[2])
    else:
        print "expr:", expr
        print "name:", name
        raise;
        return None

def EvalExpr(expr):

    ss = ""

    for e in expr:
        if isinstance(e, list):
            name = e[0]
        else:
            name = e

        if name in expTables:
            method = expTables[name]
            ss += method(e)

    return ss

def eAssign(exp):
    slots = exp[1]

    ss = ""

    if len(slots) == 0:
        # maybe assignment of matrix or struct type.

        lhs = constructExpr(exp[2])
        rhs = constructExpr(exp[3])


        if isinstance(lhs, VarRef):
            (ty, n) = GetTypeOfSymbol(lhs.orgname)
            if IsMatrixType(ty):
                
                for j in range(n):
                    for i in range(n):
                        if isinstance(rhs, Constant):
                            idx = j * n + i
                            ss += Indent() + "%s.v[%d][%d] = %s;\n" % (lhs, j, i, rhs.values[idx])
                        else:
                            ss += Indent() + "%s.v[%d][%d] = %s.v[%d][%d];\n" % (lhs, j, i, rhs, j, i)

            else:
                sym = GetSymbol(lhs.orgname)
                if sym is not None:
                    (sty, sn) = sym['type']
                    if FindStruct(sty) is not None:
                        # struct type
                        ss += Indent() + "%s = %s;\n" % (lhs, rhs)
                    else:
                        print "Invalid definition:" + sym
                        raise
                else:
                    print "Unknown or unsupported type:" + ty
                    raise
        else:
            print "Unknown assign op"
            raise

    else:

        # Don't emit code for redundant assignment 
        #
        # e.g. assign to `assignment_tmp` in global scope.
        # (assign  (xyz) (var_ref assignment_tmp)  (var_ref normalize_retval) )
        if GetScopeLevel() ==1:
            if len(exp[2]) == 2 and exp[2][0] == 'var_ref' and exp[2][1] == 'assignment_tmp':
                return "// killed redundant assign to 'assignment_tmp'\n"

        # @fixme. Supports first elem only at this time.
        slot = slots[0]

        lhs = constructExpr(exp[2])
        rhs = constructExpr(exp[3])

        # declare temp value if exist
        ss += lhs.getDeclareString() + "\n"
        ss += rhs.getDeclareString() + "\n"

        # emit intermediate expr
        ss += lhs.getIntermExprString();
        ss += rhs.getIntermExprString();

        # body
        for (i, s) in enumerate(slot):
            #print "lhs:" + lhs
            #print "rhs:" + str(rhs)
            ss += Indent() + lhs.getExprString(s, i) + " = " + rhs.getExprString(s, i) + ";\n"

    # print ss
    return ss

def eExpression(exp):
    assert 0 # todo
    print exp

def eReturn(exp):

    if len(exp) < 2:    # no argument for 'return'
        ss = Indent() + "return;\n"
    else:
        retExp = constructExpr(exp[1])

        ss = Indent() + "return " + retExp.getCExpr() + ";\n"

    return ss


def eDiscard(exp):

    ss = Indent() + "__glsl_discard(__fragment); return;\n"

    return ss


def eSwizzle(expr):
    slots = expr[1]
    args  = expr[2]
    # print expr

def eCall(expr):
    name = expr[1]

    if len(expr) < 4:
        # might be void type
        dst  = False
        args = expr[2]
    else:
        dst  = expr[2]
        args = expr[3]

    # print "dst:", dst
    # print "args:", args
    # print expr

    if dst is not False:
        # dst should be var_ref
        assert dst[0] == 'var_ref'
        dstExp = constructExpr(dst)

    isBuiltin = False
    isFuncPtrBuiltin = False

    if IsBuiltinFunction(name):
        isBuiltin = True
        prefix = "__glsl_"

        if IsBuiltinTraceFunction(name) or IsBuiltinRandomFunction(name) or IsBuiltinTextureFunction(name):
            isFuncPtrBuiltin = True

    else:
        prefix = ""

    s = ""

    # dst = func(state, a, b, c, ...)
    s += Indent() + "// { 'Call' : '" + str(args) + "' }\n"
    if dst is not False:
        s += Indent() + str(dstExp) + " = " + prefix + name + "("
    else:
        s += Indent() + prefix + name + "("

    if isFuncPtrBuiltin:
        s += "__fragment"
    elif isBuiltin:
        pass
    else:
        s += "__fragment, __state"

    for (count, arg) in enumerate(args):
        # print "arg:", arg
        exp = constructExpr(arg)
        # print "exp:", exp

        if isBuiltin:
            if isFuncPtrBuiltin:
                s += ", "
        else:
            s += ", "

        tycast = ""
        isInput = False
        isVarRef  = False

        if isinstance(exp, VarRef):
            (ty, n) = GetTypeOfSymbol(exp.orgname)

            isVarRef = True
            tycast = GetTypeCastString(exp.orgname)
            if (IsUniform(exp.orgname) or IsVaryingInput(exp.orgname)) and not IsTexture(exp.orgname):
                isInput = True

        if isInput:
            if IsVectorType(ty):
                s += "/*input:vec*/(*(" + tycast + "(" + str(exp) + ".data)))"
            else:
                s += "/*input:scalar*/" + "(*(" + tycast + "(" + str(exp) + ".data)))"
            pass
        elif isVarRef:

            if IsBuiltinVariable(exp.orgname):
                (ty, n) = GetBuiltinType(exp.orgname)

                if IsVectorType(ty):
                    tycast = "(%s%d *)" % (ty, n)
                else:
                    tycast = "(%s *)" % (ty)

                s += "/*var_ref*/" + "(*(" + tycast + str(exp) + "))"
            else:
                s += "/*var_ref*/" + str(exp)
        else:
            s += exp.getCExpr()

        if (isBuiltin and (not isFuncPtrBuiltin)) and ((len(args)-1) != count):
            s += ", "

    s += ");\n"

    return s

    # assert 0

def eIf(expr):
    # print "expr:", expr
    # print "cond:", expr[1]
    # print "then:", expr[2]
    condExpr = constructExpr(expr[1])

    # statement = [expr]
    thenStmt = expr[2]
    if len(expr[3]) > 0:    # ![]
        # print "else:", expr[3]
        elseStmt = expr[3]
    else:
        elseStmt = None
    
    # print "cond:", condExpr
    # print "then:", thenStmt
    # print "else:", elseStmt

    ss = ""

    ss = Indent() + "if ("
    ss += condExpr.getCExpr()
    ss += ") {\n"

    # then expr
    IncrementIndent()
    for e in expr[2]:
        if isinstance(e, list):
            name = e[0]
        else:
            name = e
        if name in expTables:
            method = expTables[name]
            ss += method(e)
    DecrementIndent()
    ss += "\n" + Indent() + "}"

    if elseStmt is not False:
        ss += " else {\n"
        IncrementIndent()
        for e in expr[3]:
            if isinstance(e, list):
                name = e[0]
            else:
                name = e
            if name in expTables:
                method = expTables[name]
                ss += method(e)
        DecrementIndent()
        ss += "\n" + Indent() + "}\n"
    else:
        ss += "\n"

    return ss

def eLoop(expr):
    # print "expr:", expr
    declStmt = expr[1]
    initStmt = expr[2]
    condStmt = expr[3]
    tickStmt = expr[4]
    stmt = expr[5]
    # print "stmt:", stmt 

    if len(declStmt) == 0 and len(initStmt) == 0 and len(condStmt) == 0 and len(tickStmt) == 0:
        # while loop
        ss = ""

        ss += Indent() + "while (1) {\n"
        IncrementIndent()
        ss += EvalExpr(stmt)
        DecrementIndent()
        ss += Indent() + "}\n"

    else:
        # for(i = init; i < cond; i += tick) 
        ss = ""

        assert len(declStmt) == 1
        assert len(initStmt) == 1
        assert len(condStmt) == 1
        assert len(tickStmt) == 1

        declExpr = declStmt[0]
        assert declExpr[0] == 'declare'
        initExpr = initStmt[0]
        condExpr = condStmt[0]
        tickExpr = tickStmt[0]

        decl = parseDeclare(declExpr)

        ss += Indent() + "{\n"

        ss += Indent() + "// decl = " + str(declExpr) + "\n"
        ss += EvalExpr(declStmt)

        ss += Indent() + "// init = " + str(initExpr) + "\n"
        initE = constructExpr(initExpr)
        ss += Indent() + decl['name'] + " = " + initE.getCExpr() + ";\n"

        ss += Indent() + "// cond = " + str(condExpr) + "\n"
        condE = constructExpr(condExpr)
        tickE = constructExpr(tickExpr)

        ss += Indent() + "for (; " + decl['name'] + " < " + condE.getCExpr() + "; " + decl['name'] + " += " + tickE.getCExpr() + ") {\n"
        #ss += Indent() + decl['name'] + " += " + tickE.getCExpr() + ";\n"
        IncrementIndent()

        #ss += Indent() + "if (" + decl['name'] + " >= " + condE.getCExpr() + ") {\n"
        #ss += Indent() + "  break;\n"
        #ss += Indent() + "}\n"

        for e in stmt:
            # print "e:", e

            if isinstance(e, list):
                name = e[0]
            else:
                name = e

            # Filter out tick expression
            if tickStmt:
                if name == "assign":
                    # ['assign', ['x'], ['var_ref', 'i'], ['expression', 'int', '+', ['var_ref', 'i'], ['constant', 'int', ['1']]]]
                    if e[2][1]  == decl['name']:
                        # Skip this expression
                        continue

            if name in expTables:
                method = expTables[name]
                ss += method(e)

        #ss += Indent() + "// tick = " + str(tickExpr) + "\n"
        #tickE = constructExpr(tickExpr)
        #ss += Indent() + decl['name'] + " += " + tickE.getCExpr() + ";\n"

        DecrementIndent()
        ss += Indent() + "}\n"
        ss += "}\n"

    return ss

def eBreak(expr):
    # print "break"
    return Indent() + "break;\n"

def eDeclare(exp):
    quals = exp[1]
    ty = exp[2]
    #if IsArray(ty):
    #    print "array", ty 
    offt = 3;

    (ty, n) = ParseTy(ty)

    varname = exp[offt]

    # print exp
    # print "[eDeclare] ty:", (ty, n), ", var:", varname

    isBuiltin = False
    isInOut   = False

    if ('in' in quals) or ('out' in quals) or ('uniform' in quals):
        isInOut = True
    
    if not IsBuiltinVariable(varname):
        if 'in' in quals:
            # ((ty, n), qual, name)
            AddVaryingInput(varname, ty, n, quals)
        elif 'uniform' in quals:
            # ((ty, n), qual, name)
            AddUniformInput(varname, ty, n, quals)
    else:
        isBuiltin = True

    s = ""

    # user defined global var or temporary variable needs var declaration.
    if (not isInOut) and (not isBuiltin):

        # skip redundant AST in global scope, generated from MESA's GLSL compiler.
        if varname == "assignment_tmp" and (GetScopeLevel() == 1):
            s = ""
        else:
            # Add to symbol table
            AddSymbol(varname, ty, n, quals)

            decl = VarDecl(varname, ty, n)
            s = Indent() + decl.getDeclareString();

    return s

def parseDeclare(exp):
   
    d = {}

    quals = exp[1]
    ty = exp[2]
    #if IsArray(ty):
    #    print "array", ty 

    offt = 3;

    (ty, n) = ParseTy(ty)

    varname = exp[offt]

    # print exp
    # print "ty:", (ty, n), ", var:", varname

    d['type']  = (ty, n)
    d['quals'] = quals
    d['name']  = varname 

    return d

def emitCArgs(args):

    # args = ['parameters', ...]
    if len(args) == 1:
        return ""

    s = ""
    for (n, arg) in enumerate(args[1:]):    # skip 'parameters' tag
        decl = parseDeclare(arg)

        prefix = ""
        if 'out' in decl['quals']:
            prefix = "&"    # C++ reference

        if IsVectorType(decl['type'][0]):
            s += decl['type'][0] + str(decl['type'][1]) + prefix
        elif IsMatrixType(decl['type'][0]):
            s += decl['type'][0] + str(decl['type'][1]) + prefix
        else:
            s += decl['type'][0] + prefix
        s += " "
        s += decl['name']

        if n != (len(args[1:]) - 1):
            s += ", "
            
    return s

def eFunction(exp):
    name = exp[1] 
    params = exp[2]

    # consider struct type
    p = params[1].split('__')
    if len(p) > 1:
        signature = p[0] + '__'
    else:
        signature = p[0]

    args = params[2]
    statements = params[3]

    # push symbol table stack
    gSymbolTable.append(dict())

    if IsBuiltinFunction(name):
        # No need to emit builtin function declaration.
        return ""

    isMainFunction = False

    s = ""

    if name == "main":
        # entry point function
        s += "void shader(Fragment* __fragment, FragmentState* __state)"
        isMainFunction = True

    else:
        # static function.
        argStr = emitCArgs(args)
        s += "// args = " + str(args)
        s += "\n"
        s += "static "
        s += signature  # @fixme. Support compound and vector type.  
        s += " "
        s += name
        s += " ("
        # __fragment and __state arg are required to access builtin/uniform/varying variable from user-defined function 
        s += "Fragment* __fragment, FragmentState* __state"
        if len(argStr) > 0:
            s += ", "
            s += argStr
        s += ")"

        # add function parameter to symbol table
        if len(args) > 1:   # function has arguments
            for arg in args[1:]:    # skip 'parameter' tag
                quals = arg[1]
                ty = arg[2]
                # if IsArray(ty):
                #     print "array", ty 
                offt = 3;

                (ty, n) = ParseTy(ty)

                varname = arg[offt]

                AddSymbol(varname, ty, n, quals)
        


    if len(statements) < 1:
        # Seems it has no function body.
        s += " {};\n\n";
        return s

    s += "\n{\n"

    IncrementIndent()

    # Call constructor in main function
    if isMainFunction:
        s += Indent() + "__shader_internal_constructor(__state);\n"

    for stmt in statements:
        expName = stmt[0]
        if expName in expTables:
            method = expTables[expName]
            assert method
            s += method(stmt)

    DecrementIndent()

    s += "}\n\n"

    # pop symbol table
    gSymbolTable.pop()

    return s

expTables = {
    'declare'     : eDeclare
  , 'function'    : eFunction
  , 'assign'      : eAssign
  , 'expression'  : eExpression
  , 'swiz'        : eSwizzle
  , 'call'        : eCall
  , 'if'          : eIf
  , 'loop'        : eLoop
  , 'break'       : eBreak
  , 'return'      : eReturn
  , 'discard'     : eDiscard
}

def emitEntryPoint():
    """
    Emit C entry point function definition.
    """

    s  = "extern \"C\" {\n"
    s += "void  shader(Fragment* __fragment, FragmentState* __state);\n";
    s += "int   shader_info(FragmentConfig* config);\n"
    s += "}\n"

    return s

def emitInitializer():
    """
    Emit initializer function.
    """

    qualDic = {
        'in'      : "GLSL_QUALIFIER_IN"
      , 'out'     : "GLSL_QUALIFIER_OUT"
      , 'uniform' : "GLSL_QUALIFIER_UNIFORM"
    }

    s = ""
    s += "int shader_info(FragmentConfig* config)\n"
    s += "{\n"
    IncrementIndent()

    # Uniforms
    for uniform in gUniformInputs.values():

        qualStr = "GLSL_QUALIFIER_UNIFORM"

        ty = uniform['type']
        idx = uniform['index']

        s += Indent() + "// " + str(uniform) + "\n"
        s += Indent() + "config->uniformInfos[" + str(idx) + "].type.name = \"" + ty[0] + "\";\n"
        s += Indent() + "config->uniformInfos[" + str(idx) + "].type.n    = " + str(ty[1]) + ";\n"
        s += Indent() + "config->uniformInfos[" + str(idx) + "].qualifier = " + qualStr + ";\n" # redundant?
        s += Indent() + "config->uniformInfos[" + str(idx) + "].name      = \"" + uniform['name'] + "\";\n"
        s += "\n"

    s += Indent() + "config->numUniforms = %d;\n" % (len(gUniformInputs))
    s += "\n"

    # Varyings
    for varying in gVaryingInputs.values():

        qualStr = ""
        if len(varying['quals']) == 0:
            qualStr += "GLSL_QUALIFIER_NONE"
        else:
            for (n, qual) in enumerate(varying['quals']):
                qualStr += qualDic[qual]
                if (n != (len(varying['quals']) - 1)):
                    qualStr += " | "

        ty = varying['type']
        idx = varying['index']

        s += Indent() + "// " + str(varying) + "\n"
        s += Indent() + "config->varyingInfos[" + str(idx) + "].type.name = \"" + ty[0] + "\";\n"
        s += Indent() + "config->varyingInfos[" + str(idx) + "].type.n    = " + str(ty[1]) + ";\n"
        s += Indent() + "config->varyingInfos[" + str(idx) + "].qualifier = " + qualStr + ";\n"
        s += Indent() + "config->varyingInfos[" + str(idx) + "].name      = \"" + varying['name'] + "\";\n"
        s += "\n"

    s += Indent() + "config->numVaryings = %d;\n" % (len(gVaryingInputs))
    s += "\n"

    s += Indent() + "return 0; // OK\n"

    DecrementIndent()
    s += "}\n"

    return s

def ir_to_c(input_sexp_string, opts):
    """
    Converts S-expression style IR to C/C++ code.
    """

    ir_exp = parse_sexp(input_sexp_string)

    if ('-v','') in opts:
        print "IR:" + str(ir_exp)

    s = ""

    # header
    s = "#include \"glsl_runtime.h\"\n"

    # entry point
    s += emitEntryPoint()

    #
    # pass1.1: emit struct definition
    #
    # add global scope

    s += "// --> struct definition\n"
    # print ir_exp
    for e in ir_exp:
        if isinstance(e, list):
            name = e[0]
        else:
            name = e

        if name == 'structure':
            struct = {}
            # ['structure, ['name'], ['instance'], ['N'], [fields]
            # fields = [['ty'], ['name']]
            structname = e[1][0] + "__"
            struct['name'] = structname

            members    = e[4]
            memberDefs = {}
            s += "// struct_def = " + str(e) + "\n"
            s += "typedef struct {\n"
            for member in members:
                ty   = member[0][0]
                name = member[1][0]
                memberDefs[name] = {'ty': ty, 'name': name}

                s += "    " + ty + " " + name + ";\n";
            s += "} %s;\n" % structname

            struct['members'] = memberDefs

            # Add to definition
            DefineStruct(struct)

    s += "// <-- struct definition\n"

    #
    # pass1.2: emit global variable definition
    #
    # add global scope
    gSymbolTable.append(dict())

    s += "// --> global variables\n"
    for e in ir_exp:
        if isinstance(e, list):
            name = e[0]
        else:
            name = e

        if name == 'declare':
            s += eDeclare(e)

    s += "// <-- global variables\n"


    #
    # pass2: emit global variable initializer
    #
    IncrementIndent()
    s += "static void __shader_internal_constructor(FragmentState* __state) {\n"
    for e in ir_exp:
        if isinstance(e, list):
            name = e[0]
        else:
            name = e

        if name == 'assign':
            s += eAssign(e)
    DecrementIndent()
    s += "}\n"

    #
    # pass3: body
    #
     
    for e in ir_exp:
        if isinstance(e, list):
            name = e[0]
        else:
            name = e

        # skip decl and assign in gloal scope
        if name == 'declare' or name == 'assign':
            continue

        # might be member field of struct
        if isinstance(name, list):  
            continue

        if name in expTables:
            method = expTables[name]
            s += method(e)

    # initializer
    s += emitInitializer()

    # for safety
    gSymbolTable.pop()

    return s
