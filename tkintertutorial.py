#!/usr/bin/python
# -*- codin: utf-8 -*-

"""
ZetCode Tkinter tutorial

This script shows a simple window
on the screen.

"""

from Tkinter import Tk, Frame, BOTH

class Example(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent, background="white")

        self.parent = parent

        self.initUI()

    def initUI(self):

        self.parent.title("Simple")
        self.pack(fill=BOTH, expand=1)

def main():

    root = Tk()
    root.geometry("250x150+300+300")
    app = Example(root)
    root.mainloop()

if __name__ == '__main__':
    main()
