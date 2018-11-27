#!/usr/bin/python
# -*- coding: UTF-8 -*-
'''GUI version of ReacNetGenerator
====  Usage  ====
$ reacnetgeneratorgui
'''
import os
import tkinter as tk
import tkinter.filedialog as tkfd
import webbrowser
from . import ReacNetGenerator
from ._static import _static_img


class GUI():
    '''GUI class'''

    def __init__(self):
        '''Init GUI class'''
        self._filename = ''

    def gui(self):
        '''start the GUI'''
        self._top = tk.Tk()
        self._top.title("ReacNetGenerator")
        self._top.columnconfigure(4, weight=2)

        self._filetype = tk.StringVar()
        self._filetype.set('lammpsbondfile')
        self._runhmm = tk.IntVar()
        self._runhmm.set(1)
        self._openpage = tk.IntVar()
        self._openpage.set(1)

        titleimage = tk.PhotoImage(data=_static_img['img-title'])
        self._titlelb = tk.Label(self._top, image=titleimage)
        self._titlelb.image = titleimage
        self._filenamelb = tk.Label(self._top, text="Trajectory File")
        self._filenameet = tk.Entry(self._top, width=45)
        self._openbtn = tk.Button(
            self._top, text="Select File", command=self._openfiles)
        self._filetypelb = tk.Label(self._top, text="File Type")
        self._filetyperbbond = tk.Radiobutton(
            text="LAMMPS Bond file", value="lammpsbondfile", variable=self._filetype)
        self._filetyperbtype = tk.Radiobutton(
            text="LAMMPS Dump file", value="lammpsdumpfile", variable=self._filetype)
        self._runhmmcb = tk.Checkbutton(
            text="Noise filter by HMM", variable=self._runhmm, onvalue=1, offvalue=0)
        self._openpagecb = tk.Checkbutton(
            text="Show the analysis report", variable=self._openpage, onvalue=1, offvalue=0)
        self._atomnamelb = tk.Label(self._top, text="Atomic Symbol")
        self._atomnameet = tk.Entry(self._top, width=36)
        self._atomnameet.insert(0, "C H O")
        self._runbtn = tk.Button(self._top, text="Analyze", command=self._run)

        self._titlelb.grid(row=0, column=0, columnspan=8, padx=20, pady=20)
        self._filenamelb.grid(row=1, column=0, columnspan=2, padx=5, pady=5)
        self._filenameet.grid(row=1, column=2, columnspan=5, padx=5, pady=5)
        self._openbtn.grid(row=1, column=7, columnspan=1, padx=5, pady=5)
        self._filetypelb.grid(row=2, column=0, columnspan=2, padx=5, pady=5)
        self._filetyperbbond.grid(
            row=2, column=2, columnspan=2, padx=5, pady=5)
        self._filetyperbtype.grid(
            row=2, column=4, columnspan=2, padx=5, pady=5)
        self._atomnamelb.grid(row=3, column=0, columnspan=2, padx=5, pady=5)
        self._atomnameet.grid(row=3, column=2, columnspan=3, padx=5, pady=5)
        self._runhmmcb.grid(row=2, column=6, columnspan=2, padx=5, pady=5)
        self._openpagecb.grid(row=3, column=5, columnspan=3, padx=5, pady=5)
        self._runbtn.grid(row=4, column=0, columnspan=8, padx=10, pady=10)

        self._top.mainloop()

    def _run(self):
        filename = self._filenameet.get()
        if os.path.exists(filename):
            reacnetgenerator = ReacNetGenerator(inputfilename=self._filenameet.get(
            ), inputfiletype=self._filetype.get(), atomname=self._atomnameet.get().split(), runHMM=(self._runhmm.get() == 1))
            reacnetgenerator.runanddraw()
            if self._openpage.get() == 1:
                webbrowser.open_new(os.path.abspath(
                    reacnetgenerator.resultfilename))
        else:
            tk.messagebox.showinfo("Error", "File not exsit.")

    def _openfiles(self):
        self._filename = tkfd.askopenfilename()
        if self._filename != '':
            self._filenameet.delete(0, tk.END)
            self._filenameet.insert(0, self._filename)


def gui():
    '''Open GUI version of ReacNetGenerator'''
    GUI().gui()


if __name__ == 'reacnetgenerator.gui':
    print(__doc__)
