#!/usr/bin/python

import sys
from os import path

import numpy as np
import requests
import matplotlib.pyplot as plt


from jinja2 import Template

def Json2Atoms(jsonstring):
    """Read a JSON string and return an Atoms object"""

    from ase.io.jsonio import decode
    from ase.db.row import AtomsRow

    dct = decode(jsonstring)
    row = AtomsRow(dct)

    return row.toatoms(attach_calculator=False, add_additional_information=True)

SERVER = 'https://172.23.64.223/fatman'
SERVER = 'http://127.0.0.1:5001'

COMPARE_URL    = SERVER + '/compare'
METHODS_URL    = SERVER + '/methods'
TESTS_URL      = SERVER + '/tests'
TESTRESULT_URL = SERVER + '/testresult'
RESULT_URL     = SERVER + '/results'

REPORTS_BASE_DIR = "/users/ralph/work/fatman/reports"

headersection = """
<HTML>
<HEAD>
<!--https://css-tricks.com/simple-css-row-column-highlighting/-->
<STYLE>
table {
  overflow: hidden;
}

tr:hover {
  background-color: #ffa;
}

td, th {
  position: relative;
  font-size: 11;
}
td:hover::after,
th:hover::after {
  content: "";
  position: absolute;
  background-color: #ffa;
  left: 0;
  top: -5000px;
  height: 10000px;
  width: 100%;
  z-index: -1;
}
</STYLE>
</HEAD>
<BODY><H1>Method Comparison Matrix</H1>
"""

blank_pt = {'H'  : "" , 'He' : "" , 'Li' : "" , 'Be' : "" , 'B'  : "" , 'C'  : "" , 'N'  : "" ,
            'O'  : "" , 'F'  : "" , 'Ne' : "" , 'Na' : "" , 'Mg' : "" , 'Al' : "" ,
            'Si' : "" , 'P'  : "" , 'S'  : "" , 'Cl' : "" , 'Ar' : "" , 'K'  : "" ,
            'Ca' : "" , 'Sc' : "" , 'Ti' : "" , 'V'  : "" , 'Cr' : "" , 'Mn' : "" ,
            'Fe' : "" , 'Co' : "" , 'Ni' : "" , 'Cu' : "" , 'Zn' : "" , 'Ga' : "" ,
            'Ge' : "" , 'As' : "" , 'Se' : "" , 'Br' : "" , 'Kr' : "" , 'Rb' : "" ,
            'Sr' : "" , 'Y'  : "" , 'Zr' : "" , 'Nb' : "" , 'Mo' : "" , 'Tc' : "" ,
            'Ru' : "" , 'Rh' : "" , 'Pd' : "" , 'Ag' : "" , 'Cd' : "" , 'In' : "" ,
            'Sn' : "" , 'Sb' : "" , 'Te' : "" , 'I'  : "" , 'Xe' : "" , 'Cs' : "" ,
            'Ba' : "" , 'La' : "" , 'Hf' : "" , 'Ta' : "" , 'W'  : "" , 'Re' : "" ,
            'Os' : "" , 'Ir' : "" , 'Pt' : "" , 'Au' : "" , 'Hg' : "" , 'Tl' : "" ,
            'Pb' : "" , 'Bi' : "" , 'Po' : "" , 'At' : "" , 'Rn' : "" } 

pt_template_str ="""
  <html><head><style>td{ border: 1px solid #CCCCCC; width: 30px; padding: 10px; vertical-align: center; text-align: center}</style><title>Comparing  with 1</title></head><body>
  <table style='margin-left: auto; margin-right:auto; font-size:12; table-layout: fixed; width: 1500px; padding: 2px; border-collapse: collapse'>
      <tr> <td colspan=18 style="text-align: left; font-size: 16; border: none"><h3>Comparing methods {{ m1 }} and  {{ m2 }}</td> </tr>
      <tr> 
        <td><b><b>H </b></b><br />{{ '{0:s}'.format(H) }}</td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td><b><b>He </b></b><br />{{ '{0:s}'.format(He) }}</td> </tr> 
      <tr> 
        <td><b><b>Li </b></b><br />{{ '{0:s}'.format(Li) }}</td> 
        <td><b><b>Be </b></b><br />{{ '{0:s}'.format(Be) }}</td> 
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td><b>B </b><br />{{ '{0:s}'.format(B) }}</td> 
        <td><b>C </b><br />{{ '{0:s}'.format(C) }}</td> 
        <td><b>N </b><br />{{ '{0:s}'.format(N) }}</td> 
        <td><b>O </b><br />{{ '{0:s}'.format(O) }}</td> 
        <td><b>F </b><br />{{ '{0:s}'.format(F) }}</td> 
        <td><b>Ne</b><br />{{ '{0:s}'.format(Ne) }}</td> </tr> 
      <tr> 
        <td><b>Na </b><br />{{ '{0:s}'.format(Na) }}</td> 
        <td><b>Mg </b><br />{{ '{0:s}'.format(Mg) }}</td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td> </td> 
        <td><b>Al</b><br />{{ '{0:s}'.format(Al) }}</td> 
        <td><b>Si</b><br />{{ '{0:s}'.format(Si) }}</td> 
        <td><b>P </b><br />{{ '{0:s}'.format(P) }}</td> 
        <td><b>S </b><br />{{ '{0:s}'.format(S) }}</td> 
        <td><b>Cl</b><br />{{ '{0:s}'.format(Cl) }}</td> 
        <td><b>Ar</b><br />{{ '{0:s}'.format(Ar) }}</td> </tr> 
      <tr> 
        <td><b>K </b><br />{{ '{0:s}'.format(K) }}</td> 
        <td><b>Ca</b><br />{{ '{0:s}'.format(Ca) }}</td> 
        <td><b>Sc</b><br />{{ '{0:s}'.format(Sc) }}</td> 
        <td><b>Ti</b><br />{{ '{0:s}'.format(Ti) }}</td> 
        <td><b>V </b><br />{{ '{0:s}'.format(V) }}</td> 
        <td><b>Cr</b><br />{{ '{0:s}'.format(Cr) }}</td> 
        <td><b>Mn</b><br />{{ '{0:s}'.format(Mn) }}</td> 
        <td><b>Fe</b><br />{{ '{0:s}'.format(Fe) }}</td> 
        <td><b>Co</b><br />{{ '{0:s}'.format(Co) }}</td> 
        <td><b>Ni</b><br />{{ '{0:s}'.format(Ni) }}</td> 
        <td><b>Cu</b><br />{{ '{0:s}'.format(Cu) }}</td> 
        <td><b>Zn</b><br />{{ '{0:s}'.format(Zn) }}</td> 
        <td><b>Ga</b><br />{{ '{0:s}'.format(Ga) }}</td> 
        <td><b>Ge</b><br />{{ '{0:s}'.format(Ge) }}</td> 
        <td><b>As</b><br />{{ '{0:s}'.format(As) }}</td> 
        <td><b>Se</b><br />{{ '{0:s}'.format(Se) }}</td> 
        <td><b>Br</b><br />{{ '{0:s}'.format(Br) }}</td> 
        <td><b>Kr</b><br />{{ '{0:s}'.format(Kr) }}</td> </tr> 
      <tr> 
        <td><b>Rb</b><br />{{ '{0:s}'.format(Rb) }}</td> 
        <td><b>Sr</b><br />{{ '{0:s}'.format(Sr) }}</td> 
        <td><b>Y </b><br />{{ '{0:s}'.format(Y) }}</td> 
        <td><b>Zr</b><br />{{ '{0:s}'.format(Zr) }}</td> 
        <td><b>Nb</b><br />{{ '{0:s}'.format(Nb) }}</td> 
        <td><b>Mo</b><br />{{ '{0:s}'.format(Mo) }}</td> 
        <td><b>Tc</b><br />{{ '{0:s}'.format(Tc) }}</td> 
        <td><b>Ru</b><br />{{ '{0:s}'.format(Ru) }}</td> 
        <td><b>Rh</b><br />{{ '{0:s}'.format(Rh) }}</td> 
        <td><b>Pd</b><br />{{ '{0:s}'.format(Pd) }}</td> 
        <td><b>Ag</b><br />{{ '{0:s}'.format(Ag) }}</td> 
        <td><b>Cd</b><br />{{ '{0:s}'.format(Cd) }}</td> 
        <td><b>In</b><br />{{ '{0:s}'.format(In) }}</td> 
        <td><b>Sn</b><br />{{ '{0:s}'.format(Sn) }}</td> 
        <td><b>Sb</b><br />{{ '{0:s}'.format(Sb) }}</td> 
        <td><b>Te</b><br />{{ '{0:s}'.format(Te) }}</td> 
        <td><b>I </b><br />{{ '{0:s}'.format(I) }}</td> 
        <td><b>Xe</b><br />{{ '{0:s}'.format(Xe) }}</td> </tr> 
      <tr> 
        <td><b>Cs</b><br />{{ '{0:s}'.format(Cs) }}</td> 
        <td><b>Ba</b><br />{{ '{0:s}'.format(Ba) }}</td> 
        <td><b>La</b><br />{{ '{0:s}'.format(La) }}</td> 
        <td><b>Hf</b><br />{{ '{0:s}'.format(Hf) }}</td> 
        <td><b>Ta</b><br />{{ '{0:s}'.format(Ta) }}</td> 
        <td><b>W </b><br />{{ '{0:s}'.format(W) }}</td> 
        <td><b>Re</b><br />{{ '{0:s}'.format(Re) }}</td> 
        <td><b>Os</b><br />{{ '{0:s}'.format(Os) }}</td> 
        <td><b>Ir</b><br />{{ '{0:s}'.format(Ir) }}</td> 
        <td><b>Pt</b><br />{{ '{0:s}'.format(Pt) }}</td> 
        <td><b>Au</b><br />{{ '{0:s}'.format(Au) }}</td> 
        <td><b>Hg</b><br />{{ '{0:s}'.format(Hg) }}</td> 
        <td><b>Tl</b><br />{{ '{0:s}'.format(Tl) }}</td> 
        <td><b>Pb</b><br />{{ '{0:s}'.format(Pb) }}</td> 
        <td><b>Bi</b><br />{{ '{0:s}'.format(Bi) }}</td> 
        <td><b>Po</b><br />{{ '{0:s}'.format(Po) }}</td> 
        <td><b>At</b><br />{{ '{0:s}'.format(At) }}</td> 
        <td><b>Rn</b><br />{{ '{0:s}'.format(Rn) }}</td> </tr> 
      <tr> <td colspan=18 style="text-align: left; font-size: 14; border: none"><a href="{{ backlink }}">List View</a>&nbsp;&nbsp;&nbsp;&nbsp;<a href="index.html">Home</a>  </td></tr>
  </table>
    </body></html>
"""
pt_template = Template(pt_template_str)


def compare(method1, method2, test=None):
    r= requests.get(COMPARE_URL, data={'method1':method1, 'method2':method2, 'test':test}, verify=False)
    return r.json()

def eos(V0,B0, B1, E0=0.):
    B0 = B0 * 1e9 / 1.602176565e-19 / 1e30
    rng = np.linspace(0.93*V0, 1.07*V0, 40)
    E = [ E0 + 9./16. * V0 * B0 *  ( ((V0/v)**(2./3.) -1)**3 * B1 +
                                     ((V0/v)**(2./3.) -1)**2 * (6-4*(V0/v)**(2./3.)) ) for v in rng]
    return rng, np.array(E)

def tresult(method, test):
    r= requests.get(TESTRESULT_URL, data={'method':method, 'test':test}, verify=False)
    return r.json()
    
def result(method, test=None):
    r= requests.get(RESULT_URL, data={'method':method, 'test':test}, verify=False)
    return r.json()
    
def energies_for_test(method, test):
    V = []
    E = []
    r = result(method=method,test=test)
    for x in r:
        natom = len(Json2Atoms(x["task"]['structure']['ase_structure']).get_masses())
        V.append(Json2Atoms(x["task"]['structure']['ase_structure']).get_volume()/natom)
        E.append(x['energy']/natom)
        
    return np.array(V),np.array(E)

def full_plot_data(method, test=None):
    a = tresult(method, test)
    V0 = a[test]['V']
    B0 = a[test]['B0']
    B1 = a[test]['B1']
    
        
    xfit,yfit = eos(V0,B0,B1)
    
    x2,y2 = energies_for_test(method=method, test=test)
    if 'E0' in a[test].keys():
        E0 = a[test]['E0']
        y2 -= E0
    
    return [xfit,yfit, x2,y2]

def full_compare(method1, method2, writefile = False):
    r = requests.get(COMPARE_URL, params={'method1':method1, 'method2':method2}, verify = False)
    f = r.json()
    tests = f['test'].keys()
    
    for t in tests:
        c = full_plot_data(method1,t)
        r = full_plot_data(method2,t)
        if len(c[2])==0 and len(r[2]) == 0 :
            volpoints = []
            epoints = []
        elif len(c[2])>0:
            volpoints = c[2]
            epoints = c[3]
        else:
            volpoints = r[2]
            epoints = r[3]
        
        myfig = plt.figure(figsize=(15,4))
        ax = plt.subplot2grid((2,2),(0, 0),rowspan=2)
        ax.plot(c[0],c[1],'r', c[2],c[3], 'rs', label=['calc','fitted'])
        ax.plot(r[0],r[1],'b', r[2],r[3], 'bs', label=['calc','fitted'])
        ax.annotate("Method {:}".format(method1), xy= (c[0][3], c[1][3]), xytext = (-5,-30), textcoords = "offset points", fontsize=10, arrowprops=dict(facecolor='black', shrink=0.05, headwidth=3, width=1), horizontalalignment='right')
        ax.annotate("Method {:}".format(method2), xy= (r[0][-3], r[1][-3]), xytext = (5,-30), textcoords = "offset points", fontsize=10, arrowprops=dict(facecolor='black', shrink=0.05, headwidth=3, width=1), horizontalalignment='left')
        diffc = np.interp(c[0],r[0],r[1])
        #ax.fill_between(c[0],c[1],diffc, facecolor="#999999", alpha=0.7, color="#999999")
        ax.fill_between(c[0], c[1], diffc, where=c[1] >= diffc, facecolor='#111111',alpha=0.6, color="#111111")
        ax.fill_between(c[0], c[1], diffc, where=c[1] <= diffc, facecolor='#AAAAAA',alpha=0.6, color="#AAAAAA")
        
        ax.set_title(t)
        
        ax = plt.subplot2grid((2,2),(0, 1))
        ax.plot(c[0],c[1]-r[1],'r')
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_title("Residual")
        if not len(volpoints)==0:
            ax.plot(volpoints,[0 for x in volpoints], 'ks', markersize=2.5)
            xtmp = np.linspace(volpoints[2]*0.94, volpoints[2]*1.06, 40)
            ytmp = np.interp(xtmp, c[0], c[1]-r[1])
            #ax.axvline(xtmp[0], color='k', linewidth=0.5)
            #ax.axvline(xtmp[-1], color='k', linewidth=0.5)
            ax.fill_between(xtmp, [0. for x in xtmp], ytmp, where=0 >= ytmp, facecolor='#111111',alpha=0.6, color="#AAAAAA")
            ax.fill_between(xtmp, [0. for x in xtmp], ytmp, where=0 <= ytmp, facecolor='#AAAAAA',alpha=0.6, color="#AAAAAA")

        ax = plt.subplot2grid((2,2),(1, 1))
        ax.text(0.4,0.5,"$\Delta V_0 = {:>8.5f}\quad   ({:4.3f}\%)$".format(f['test'][t][3]-f['test'][t][0],(f['test'][t][3]-f['test'][t][0])/f['test'][t][3]*100) ,fontsize=16)
        ax.text(0.4,0.3,"$\Delta B_0 = {:>8.5f}\quad   ({:4.3f}\%)$".format(f['test'][t][4]-f['test'][t][1],(f['test'][t][4]-f['test'][t][1])/f['test'][t][4]*100) ,fontsize=16)
        ax.text(0.4,0.1,"$\Delta B_1 = {:>8.5f}\quad   ({:4.3f}\%)$".format(f['test'][t][5]-f['test'][t][2],(f['test'][t][5]-f['test'][t][2])/f['test'][t][5]*100) ,fontsize=16)
        ax.text(0,0.3,"$\mathbf{{\Delta = {:8.3f}}}\,\mathrm{{meV}}$".format(f['test'][t][6]),fontsize=18)
        
        ax.axis("off")
        if writefile: 
            if "deltatest_" in t:
                element = t.replace("deltatest_","")
            
            picture = "img/{:04d}_{:04d}_{:s}.png".format(method1,method2,element)

            plt.savefig(path.join(REPORTS_BASE_DIR, "html", picture), dpi=200)
        plt.close(myfig)


class deltaReport():
    """ Base Class for reporting Deltatest results """
    def __init__(self, outfile):
        self.outfile = outfile
        self.output = ""
        self.nlines = 0

    def set_report_header(self, code="abinit", subtitle="small text", features=[]):
        self.output = self.get_file_intro() 
        self.output += self.get_method_description(code, subtitle, features)


    def set_table_footer(self, delta_avg=-1., delta_std=-1.):
        self.output += self.get_table_outro(delta_avg, delta_std) 

    def set_report_footer(self):
        self.output += self.get_file_outro() 


    def write(self):
        with open(self.outfile,'w') as of:
            of.write(self.output)



class HTMLReport(deltaReport):
    """Output as fancy HTML table"""
    def __init__(self, outfile):
        deltaReport.__init__(self, outfile)
        self.added_headers = False
        self.added_link_headers = False
        self.columns = 0

    def get_file_intro(self):
        text = r"""<HTML><HEAD>
                         <TITLE>Deltatest Report</TITLE>
                        <script type="text/javascript" src="jquery-latest.js"></script> 
                        <script type="text/javascript" src="jquery.tablesorter.min.js"></script> 
                        <link rel="stylesheet" type="text/css" href="blue/style.css">
                        <SCRIPT type="text/javascript" >
                         $(document).ready(function() 
                            { 
                                $("table").tablesorter( {sortList: [[0,0]]} ); 
                            } 
                        );  
                        </SCRIPT>
                        </HEAD><BODY>
                """
        return text


    def get_method_description(self, code, subtitle, features):
        ncols = len(features)/2+1
        line1 = ""
        line2 = ""

        for f,v in features[::2]:
            line1 += r"<TD> <b>{:s}</b>: {:s} </TD>".format(f,v)
        for f,v in features[1::2]:
            line2 += r"<TD> <b>{:s}</b>: {:s} </TD>".format(f,v)

        text  = r"""<TABLE width="1100px" frame="box" align="center">"""
        text += r"""<TR><TD><font size=4> {:s}</font></TD> {:s}</TR>""".format(code,line1)
        text += r"""<TR><TD>{:s}</TD> {:s}</TR>""".format(subtitle, line2)


        text += r"""</TABLE>"""
        return text

    def add_link(self, url, disp=None):
        if not self.added_link_headers:
            self.add_link_headers()

        if disp is None:
            disp=url
        text = """
                   &nbsp;&nbsp;<A href="{0:s}">{1:s}</A>&nbsp;&nbsp;\n
               """.format(url,disp)
        self.output += text

    def add_link_headers(self):
        text = """
                   <TABLE width="1100px" align="center"><TR><TD>Subgroups:
               """
        self.output += text
        self.added_link_headers = True


    def get_file_outro(self):
        text = r""" </BODY>
                    </HTML>
                    """
        if self.added_link_headers:
            text = "</TD></TR></TABLE>" + text
        return text

    def get_table_outro(self, delta_avg=None, delta_std=None):
        if delta_avg is not None:
            if delta_std is not None:
                avgtext = "<TR style=\"border-top:1px solid;border-bottom:1px solid;\"><TD colspan=\"{:d}\" align=\"right\">Average Delta (n={:2d}):</TD><TD align=\"right\">{:6.3f}</TD><TD align=\"right\">&plusmn; {:6.3f}</TD></TR>".format(self.columns-2, self.nlines, delta_avg,delta_std)   
            else:
                avgtext = "<TR><TD colspan=\"{:d}\" align=\"right\">Average Delta (n={:2d}):</TD><TD> {:6.3f}</TD></TR>".format(self.columns-1, self.nlines, delta_avg)   
        else:
            avgtext = ""

        text = r""" </TBODY>{:s}
                    </TABLE>
                """.format(avgtext)
        return text

    def add_line(self, dataline):
        if self.added_headers == False:
            self.output += self.get_tableheaders(dataline)
            self.added_headers = True
            self.columns = len(dataline)

        self.nlines+=1
        text = "<TR>"
        data = [x[1] for x in dataline]

        for d in data:
            if isinstance(d, float):
                text += "<TD align=\"right\">{:6.3f}</TD>".format(d)
            elif isinstance(d, tuple):
                text += "<TD align=\"right\">{:>2d}, {:>2d}, {:>2d}</TD>".format(d[0], d[1], d[2])
            elif isinstance(d, str):
                text += "<TD align=\"right\">{:s}</TD>".format(d)
            else:
                text += "<TD align=\"right\">{:s}</TD>".format(str(d))

        text = text[:-1] + "</TR>\n"
        self.output  += text

    def get_tableheaders(self, dataline):
        text = "\n\n <br/> <TABLE width=\"1100px\" align=\"center\" id=\"resultsTable\" class=\"tablesorter\" style=\"border-collapse:collapse;\"> \n<THEAD><TR style=\"border-top:1px solid;border-bottom:1px solid;\">"
        keys = [x[0] for x in dataline]

        for k in keys:
            text += "<TH align=\"right\"><b>{:s}</b></TH>".format(k)

        text = text + "</TR></THEAD><TBODY> \n"
        return text


def create_html_comparison(idlist = None):
    """make a bunch of comparative html outputs"""
    elements =  { "H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7, "O":8, "F":9, "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15, "S":16, "Cl":17, "Ar":18, "K":19, "Ca":20, "Sc":21, "Ti":22, "V":23, "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28, "Cu":29, "Zn":30, "Ga":31, "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36, "Rb":37, "Sr":38, "Y":39, "Zr":40, "Nb":41, "Mo":42, "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47, "Cd":48, "In":49, "Sn":50, "Sb":51,  "Te":52, "I":53, "Xe":54, "Cs":55, "Ba":56, "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77, "Pt":78, "Au":79, "Hg":80, "Tl":81, "Pb":82, "Bi":83,  "Po":84, "Rn":86 }

    req = requests.get(METHODS_URL, verify = False)
    req.raise_for_status()
    method_list = sorted(req.json(), key = lambda x: x['id'])

    of = open(path.join(REPORTS_BASE_DIR, "html/index.html"), "w")
    of.write(headersection)
    of.write("<TABLE><TR><TD style=\"width:40px\"></TD>")
    of.write("".join(["<TD style=\"width:40px\"><span title=\"ID: {id}, code: {code}, pseudopotential: {pseudopotential}, basis set: {basis_set}\">{id}</span></TD>".format(**m) for m  in method_list]))
    of.write("</TR>\n")

    for method_row in method_list:
        desc_1 = "ID: {id}, code: {code}, pseudopotential: {pseudopotential}, basis set: {basis_set}".format(**method_row)
        of.write("<TR style=\"hover {{background: yellow;}}\"><TD><span title=\"{desc}\">{id}</span></TD>".format(desc=desc_1, **method_row))

        m_id_1 = method_row['id']
        for method_col in method_list:
            m_id_2 = method_col['id']
            desc_2 = "ID: {id}, code: {code}, pseudopotential: {pseudopotential}, basis set: {basis_set}".format(**method_col)
            if m_id_2 <= m_id_1:
                of.write("<TD></TD>")
                continue
            print('{} {}'.format(m_id_1, m_id_2))

            req = requests.get(COMPARE_URL, params = {"method1": m_id_1, "method2": m_id_2}, verify=False)

            req.raise_for_status()
            a = req.json()
            val = a["summary"]["avg"]
            if val > 99:
                of.write("<TD><a href={:} title=\"avg: {:}\n std: {:}\n n: {:}\">&gt;99</a></TD>".format("{:04d}-{:04d}.html".format(m_id_1, m_id_2),val,a["summary"]["stdev"],a["summary"]["N"]))
            elif str(val)=='nan':
                of.write("<TD></TD>")
                continue
            else:
                of.write("<TD><a href={:} title=\"avg: {:}\n std: {:}\n n: {:}\">{:3.2f}</a></TD>".format("{:04d}-{:04d}.html".format(m_id_1, m_id_2),val,a["summary"]["stdev"],a["summary"]["N"], val ))

            #skipping non-requested reports.
            if idlist is not None and not (m_id_1 in idlist or m_id_2 in idlist):
                continue

            full_compare(m_id_1, m_id_2, writefile=True)

            pt_results = blank_pt.copy()
            pt_results['m1'] = m_id_1
            pt_results['m2'] = m_id_2
            pt_results['backlink'] = "{:04d}-{:04d}.html".format(m_id_1, m_id_2)

            detailreport = HTMLReport(path.join(REPORTS_BASE_DIR, "html", "{:04d}-{:04d}.html".format(m_id_1, m_id_2)))
            detailreport.set_report_header(code=desc_1, subtitle=desc_2, features=[('View', '<a href="{:}">PT</a>'.format("{:04d}-{:04d}-pt.html".format(m_id_1, m_id_2)))])
            #detailreport.
            for t, line in a["test"].items():
                if "deltatest_" in t:
                    element = t.replace("deltatest_","")
                    pt_results[element] = '{:4.3f}'.format(line[6])
                
                picture = "img/{:04d}_{:04d}_{:s}.png".format(m_id_1, m_id_2, element)

                dataline = [("z", str(elements[element])),
                            ("Element", "<a href='{:}'>{:}</a>".format(picture, element)),
                            ("V<sub>0</sub>", line[0]), 
                            ("B<sub>0</sub>", line[1]),  
                            ("B<sub>1</sub>", line[2]), 
                            ("V<sub>0,r</sub>", line[3]), 
                            ("B<sub>0,r</sub>", line[4]),  
                            ("B<sub>1,r</sub>", line[5]), 
                            ("&Delta;",         line[6]) ]
                detailreport.add_line(dataline)

            detailreport.set_table_footer(delta_avg = a["summary"]["avg"], delta_std=a["summary"]["stdev"])
            detailreport.write()

            with open(path.join(REPORTS_BASE_DIR, "html", "{:04d}-{:04d}-pt.html".format(m_id_1, m_id_2)),'w') as pt_outfile:
                pt_outfile.write(pt_template.render(pt_results))

        of.write("</TR>")


    of.write("</TABLE>")
    of.write("</BODY></HTML>")
    of.close()









    req = requests.get(TESTS_URL, verify = False)
    req.raise_for_status()
    test_list = sorted(req.json(), key = lambda x:x[0])

    of = open(path.join(REPORTS_BASE_DIR, "html-by-test/index.html"), "w")
    of.write(headersection)
    of.write("<TABLE><TR><TD style=\"width:40px\"></TD>")
    of.write("".join(["<TD style=\"width:40px\"><span title=\"ID: {id}, code: {code}, pseudopotential: {pseudopotential}, basis set: {basis_set}\">{id}</span></TD>".format(**m) for m  in method_list]))
    of.write("</TR>\n")

    for t_id_1, desc_1 in test_list:
        if not 'deltatest' in desc_1: continue
        of.write("<TR style=\"hover {{background: yellow;}}\"><TD>{:}</TD>".format(desc_1))

        detailreport = HTMLReport(path.join(REPORTS_BASE_DIR, "html-by-test", "{:}.html".format(desc_1)))
        detailreport.set_report_header(code=desc_1, subtitle="Reference V0 = {:}, B0 = {:}, B1 = {:}", features=[])

        for method in method_list:
            m_id_2 = method['id']
           #if m_id_2<=m_id_1 :
           #    of.write("<TD></TD>")
           #    continue
            print('{} {}'.format(desc_1, m_id_2))
            req = requests.get(COMPARE_URL, params = {"method1": m_id_2, "method2": 3, "test": desc_1} , verify = False)
            req.raise_for_status()
            a = req.json()
            val = a["summary"]["avg"]
            if val>99:
                of.write("<TD><a href={:}.html>&gt;99</a></TD>".format(desc_1))
            elif str(val)=='nan':
                of.write("<TD></TD>")
                continue
            else:
                of.write("<TD><a href={:}.html>{:3.2f}</a></TD>".format(desc_1, val))

            #detailreport.
            for t, line in a["test"].items():
                if "deltatest_" in t:
                    element = t.replace("deltatest_","")
                
                no1 = 3 if m_id_2 > 3 else m_id_2
                no2 = m_id_2 if m_id_2 >= 3 else 3

                picture = "img/{:04d}_{:04d}_{:s}.png".format(no1, no2, element)

                dataline = [#("z", str(elements[element])),
                            ("ID", method['id']),
                           #("Element", element), 
                            ("Method", "<a href='{picture}'>ID: {id}, code: {code}, pseudopotential: {pseudopotential}, basis set: {basis_set}</a>".format(picture=picture, **method)),
                            ("V<sub>0</sub>", line[0]), 
                            ("B<sub>0</sub>", line[1]),  
                            ("B<sub>1</sub>", line[2]), 
                           #("V<sub>0,r</sub>", line[3]), 
                           #("B<sub>0,r</sub>", line[4]),  
                           #("B<sub>1,r</sub>", line[5]), 
                            ("&Delta;",         line[6]),
                            ]
                detailreport.add_line(dataline)

        detailreport.set_table_footer(delta_avg = 0., delta_std=0.)
        detailreport.write()
        of.write("</TR>")


    of.write("</TABLE>")
    of.write("</BODY></HTML>")
    of.close()

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 0:
        idlist = [int(x) for x in args]
        create_html_comparison(idlist)
    else:
        create_html_comparison()


