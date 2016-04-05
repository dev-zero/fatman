#!/usr/bin/python
import os
import numpy as np
import requests

SERVER = 'https://172.23.64.223/fatman'
SERVER = 'http://127.0.0.1:5001'
COMPARE_URL = SERVER + '/compare'
METHODS_URL = SERVER + '/methods'


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
                text += "<TD align=\"right\">{:s}</TD>".format(d)

        text = text[:-1] + "</TR>\n"
        self.output  += text

    def get_tableheaders(self, dataline):
        text = "\n\n <br/> <TABLE width=\"1100px\" align=\"center\" id=\"resultsTable\" class=\"tablesorter\" style=\"border-collapse:collapse;\"> \n<THEAD><TR style=\"border-top:1px solid;border-bottom:1px solid;\">"
        keys = [x[0] for x in dataline]

        for k in keys:
            text += "<TH align=\"right\"><b>{:s}</b></TH>".format(k)

        text = text + "</TR></THEAD><TBODY> \n"
        return text


def create_html_comparison():
    """make a bunch of comparative html outputs"""

    req = requests.get(METHODS_URL, verify = False)
    req.raise_for_status()
    method_list = sorted(req.json(), key = lambda x:x[0])

    of = open("/users/ralph/work/fatman/reports/html/index.html","w")
    of.write("<HTML><BODY><H1>REPORT OVERVIEW</H1>\n")
    of.write("<TABLE><TR><TD style=\"width:40px\"></TD>")
    of.write("".join(["<TD style=\"width:40px\"><span title=\"{:}\">{:}</span></TD>".format(x[1],x[0]) for x  in method_list]))
    of.write("</TR>")

    for m_id_1, desc_1 in method_list:
        of.write("<TR><TD><span title=\"{:}\">{:}</span></TD>".format(desc_1,m_id_1))

        for m_id_2, desc_2 in method_list:
            req = requests.get(COMPARE_URL, params = {"method1": m_id_1, "method2": m_id_2} , verify = False)

            req.raise_for_status()
            a = req.json()
            val = a["summary"]["avg"]
            if val>99:
                of.write("<TD>&gt;99</TD>")
            elif val==np.nan:
                of.write("<TD></TD>")
            else:
                of.write("<TD><span title=\"avg: {:}\n std: {:}\n n: {:}\">{:3.2f}</span></TD>".format(val,a["summary"]["stdev"],a["summary"]["N"], val ))
            print a["methods"], a["summary"]["avg"]
            #quit()
        of.write("</TR>")


    of.write("</TABLE>")
    of.write("</BODY></HTML>")
    of.close()

if __name__ == "__main__":
    create_html_comparison()


