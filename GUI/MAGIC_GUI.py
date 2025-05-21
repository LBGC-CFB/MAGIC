import os
import sys
import tkinter as tk
import tkinter.font
from tkinter import ttk
from tkinter import messagebox
from tkinter import Label, Entry, Button, StringVar, IntVar, PhotoImage, SUNKEN
from tkinter.filedialog import askdirectory, askopenfilename
from tkinter.scrolledtext import ScrolledText



class MAGIC_GUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)
        # define path of ico and ppm
        base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
        icon_path = os.path.join(base_path, "MAGIC.ico")
        ppm_path = os.path.join(base_path, "MAGIC.ppm")

        # define main window
        self.title("MAGIC, create reference files")
        self.iconbitmap(icon_path)
        self.geometry('1500x900')
        self.resizable(True, True)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        # define canvas and frame (used for scrollbar)
        self.frame = ttk.Frame(self)
        self.frame.grid(row=0, column=0, sticky="nsew")

        self.canvas = tk.Canvas(self.frame)
        self.scrollbarv = ttk.Scrollbar(self.frame, orient="vertical", command=self.canvas.yview)
        self.scrollbarh = ttk.Scrollbar(self.frame, orient="horizontal", command=self.canvas.xview)

        self.scrollable_frame = ttk.Frame(self.canvas)
        self.scrollable_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbarv.set)
        self.canvas.configure(xscrollcommand=self.scrollbarh.set)

        self.scrollbarv.grid(row=0, column=1, sticky="ns")
        self.scrollbarh.grid(row=1, column=0, sticky="ew")
        self.canvas.grid(row=0, column=0, sticky="nsew")
    
        self.frame.grid(row=0, column=0, sticky="nsew")
        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(0, weight=1)
        self.canvas.bind_all("<MouseWheel>", self._on_mousewheel)

        # define logo and attach to previously created frame
        self.logo = PhotoImage(file=ppm_path)
        self.logo_label = Label(self.scrollable_frame, image=self.logo)
        self.logo_label.grid(row=0, column=1, columnspan=6, pady=20) #, sticky="nsew")

        # define tkinter variables
        self.outdir = StringVar()
        self.seq = StringVar()
        self.UTR5p = StringVar()
        self.UTR3p = StringVar()
        self.gene_name = StringVar()
        self.transcript_name = StringVar()
        self.construction_name = StringVar()
        self.exon_number = IntVar()
        self.exs_list = []

        # define fonts
        self.ask_font = tkinter.font.Font(family="Helvetica", size="10", weight='bold')
        self.answer_font = tkinter.font.Font(family="Helvetica", size="10")
        self.entry_font = tkinter.font.Font(family="Helvetica", size="10")
        self.button_font = tkinter.font.Font(family="Helvetica", size="10")
        self.button_exit_font = tkinter.font.Font(family="Helvetica", size="17", weight='bold')

        # build gui
        self.add_widgets()


    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        #self.canvas.xview_scroll(int(-1 * (event.delta / 120)), "units")


    def add_widgets(self):
        # widget1: ask output directory
        self.ask_dir = Label(self.scrollable_frame, text='Select ouput directory', justify="left", font=self.ask_font)
        self.ask_dir.grid(row=5, column=0, padx=5, pady=10, ipadx=3, ipady=3, sticky="w")
        self.ask_btn_dir = Button(self.scrollable_frame, text ='Choose directory', command = lambda:self.choose_outdir(), font=self.button_font)
        self.ask_btn_dir.grid(row=5, column=1, ipadx=1, ipady=1)

        # widget2: ask fasta file
        self.ask_input = Label(self.scrollable_frame, text='Select input file', justify="left", font=self.ask_font)
        self.ask_input.grid(row=6, column=0, padx=5, pady=10, ipadx=3, ipady=3, sticky="w")
        self.ask_btn_input = Button(self.scrollable_frame, text ='Choose input file', command = lambda:self.choose_fasta_file(), font=self.button_font)
        self.ask_btn_input.grid(row=6, column=1, ipadx=1, ipady=1)

        # widget3: ask exon-specific vector sequences 5' and 3'
        self.UTR_seqs = []
        self.ask_seq_ex_5p = Label(self.scrollable_frame, text='Enter: \nconstitutive exon vector \n5 prime sequence', justify="left", font=self.ask_font)
        self.ask_seq_ex_5p.grid(row=7, column=0, padx=5, pady=10, sticky="w")
        self.entry_seq_ex_5p = ScrolledText(self.scrollable_frame, font=self.answer_font, height=1, width=35)
        self.entry_seq_ex_5p.grid(row=7, column=1, padx=5, pady=10, ipadx=3, ipady=3)#, columnspan=3)

        self.ask_seq_ex_3p = Label(self.scrollable_frame, text='constitutive exon vector \n3 prime sequence', justify="left", font=self.ask_font)
        self.ask_seq_ex_3p.grid(row=7, column=2, padx=5, pady=10, sticky="w")
        self.entry_seq_ex_3p = ScrolledText(self.scrollable_frame, font=self.answer_font, height=1, width=35)
        self.entry_seq_ex_3p.grid(row=7, column=3, padx=5, pady=10, ipadx=3, ipady=3)#, columnspan=3)

        self.ask_btn_seq_ex_5p3p = Button(self.scrollable_frame, text='Validate ✅', command = lambda:self.get_UTR_informations(), bg="green", fg="white", font=self.ask_font)
        self.ask_btn_seq_ex_5p3p.grid(row=7, column=4, ipadx=1, ipady=1, sticky="e")

        # widget4, 5 & 6: ask gene, transcript and construction name
        self.ask_gene = Label(self.scrollable_frame, text=' Enter: Gene name', justify="left", font=self.ask_font)
        self.ask_gene.grid(row=8, column=0, padx=5, pady=10, sticky="w")
        self.entry_gene = Entry(self.scrollable_frame, text='Enter: Gene name', font=self.answer_font)
        self.entry_gene.grid(row=8, column=1, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_tr = Label(self.scrollable_frame, text='Transcript name', justify="left", font=self.ask_font)
        self.ask_tr.grid(row=8, column=2, padx=5, pady=10, ipadx=3, ipady=3)
        self.entry_tr = Entry(self.scrollable_frame, text='Transcript name', font=self.answer_font)
        self.entry_tr.grid(row=8, column=3, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_const = Label(self.scrollable_frame, text='Construction name', justify="left", font=self.ask_font)
        self.ask_const.grid(row=8, column=4, padx=5, pady=10, ipadx=3, ipady=3)
        self.entry_const = Entry(self.scrollable_frame, text='Construction name', font=self.answer_font)
        self.entry_const.grid(row=8, column=5, padx=5, pady=10, ipadx=3, ipady=3)

        self.ask_btn_gene_tr_const = Button(self.scrollable_frame, text='Validate ✅', command = lambda:self.get_construction_informations(), bg="green", fg="white", font=self.ask_font)
        self.ask_btn_gene_tr_const.grid(row=8, column=6, ipadx=1, ipady=1, sticky="e")

        # widget 7: ask exon(s) number(s)
        self.ask_ex = Label(self.scrollable_frame, text='Select number of exons', font=self.ask_font)
        self.ask_ex.grid(row=9, column=0, padx=5, pady=10, ipadx=3, ipady=3)
        self.ask_btn_ex = Button(self.scrollable_frame, text ='Select number of exons', command = lambda:self.ask_exon(), font=self.button_font)
        self.ask_btn_ex.grid(row=9, column=1, ipadx=1, ipady=1)

        # final exit button to generate files
        self.exit_button = Button(self.scrollable_frame, text="★⋆⭒˚.⋆  Generate files  ⋆⭒˚.⋆★", command=self.destroy, font=self.button_exit_font, bg="#EC0C6E", fg="#F8D7E5").grid(row=1000, column=2, pady=20, columnspan=3)


    def choose_outdir(self):
        self.selected_dir = askdirectory()
        if self.selected_dir is not None:
            tk.Label(self.scrollable_frame, text=f"{self.selected_dir}", font=self.answer_font, fg="#8D064B").grid(row=5, column=2, padx=5, pady=10, columnspan=4, sticky="w")
            self.outdir.set(self.selected_dir)
            self.ask_btn_dir.config(relief=SUNKEN, fg="gray")


    def choose_fasta_file(self):
        self.file_path = askopenfilename(filetypes=([("all_files","*.fa"), ("all_files","*.fasta"), ("all_files","*.gb")]))
        if self.file_path is not None:
            tk.Label(self.scrollable_frame, text=f"{self.file_path}", font=self.answer_font, fg="#8D064B").grid(row=6, column=2, padx=5, pady=10, columnspan=4, sticky="w")
            if self.file_path.split(".")[-1] == "gb":
                seq_list = []
                with open (self.file_path, "r") as filin:
                    for line in filin:
                        line = line.split()
                        if line[0].isdigit():
                            for ele in line[1:]:
                                seq_list.append(ele)
                self.seq.set("".join(seq_list))
            elif self.file_path.split(".")[-1] == "fasta":
                with open (self.file_path, "r") as filin:
                    for line in filin:
                        if not line.startswith(">"):
                            self.seq.set(line.strip())
            tk.Label(self.scrollable_frame, text=f"sequence length: {len(self.seq.get())} bp", font=self.answer_font, fg="#8D064B").grid(row=6, column=5, padx=5, pady=10)
            self.ask_btn_input.config(relief=SUNKEN, fg="gray")


    def get_UTR_informations(self):
        self.UTR5p.set(self.entry_seq_ex_5p.get("1.0", tk.END).strip())
        self.UTR3p.set(self.entry_seq_ex_3p.get("1.0", tk.END).strip())
        self.ask_btn_seq_ex_5p3p.config(relief=SUNKEN, bg="white", fg="green")


    def get_construction_informations(self):    
        self.gene_name.set(self.entry_gene.get())
        self.transcript_name.set(self.entry_tr.get())
        self.construction_name.set(self.entry_const.get())
        self.ask_btn_gene_tr_const.config(relief=SUNKEN, bg="white", fg="green")


    def ask_exon(self):
        self.exons = tk.simpledialog.askinteger("Number of exons?", "Number of exons?", parent = self)
        self.ask_btn_ex.config(relief=SUNKEN, state="disabled")
        if self.exons:
            tk.Label(self.scrollable_frame, text=f"Exon(s): {self.exons}", font=self.answer_font, fg="#8D064B").grid(row=9, column=2, padx=5, pady=5)
            self.ex_infos = []
            for i in range(self.exons):
                # widget 1: ask exon number
                self.ask_exnb = Label(self.scrollable_frame, text=f'Enter: Exon {i+1} number?', font=self.ask_font)
                self.ask_exnb.grid(row=10+i+1, column=0, padx=5, pady=5, ipadx=3, ipady=3, sticky="w")
                self.entry_exnb = Entry(self.scrollable_frame, font=self.answer_font)
                self.entry_exnb.grid(row=10+i+1, column=1, padx=5, pady=5, ipadx=3, ipady=3)

                # widget 2: ask exon sequence
                self.ask_exseq = Label(self.scrollable_frame, text=f'Exon {i+1} sequence?', font=self.ask_font)
                self.ask_exseq.grid(row=10+i+1, column=2, padx=5, pady=5, ipadx=3, ipady=3, sticky="w")
                self.stext_exseq = ScrolledText(self.scrollable_frame, font=self.answer_font, height=2, width=70)
                self.stext_exseq.grid(row=10+i+1, column=3,  padx=5, pady=5, ipadx=3, ipady=3, columnspan=4)
                self.ex_infos.append([self.entry_exnb, self.stext_exseq])

            self.ask_ex_num_seq = Button(self.scrollable_frame, text="Validate ✅", command=lambda:self.keep_exons_info(lst=self.exs_list), bg="green", fg="white", font=self.ask_font)
            self.ask_ex_num_seq.grid(row=10+i+1, column=8, padx=5, pady=5, ipadx=1, ipady=1)


    def keep_exons_info(self, lst):
        lst.clear()
        for ex_inf in self.ex_infos:
            tmp_lst = []
            tmp_lst.append(ex_inf[0].get())
            tmp_lst.append(ex_inf[1].get("1.0","end-1c"))
            lst.append(tmp_lst)
        self.ask_ex_num_seq.config(relief=SUNKEN, bg="white", fg="green")

        

class generate_files():
    def __init__(self):
        # retrieve value from MAGIC GUI
        self.gene_name = magic.gene_name.get()
        self.outdir = magic.outdir.get()
        self.seq = magic.seq.get()
        self.UTR5p = magic.UTR5p.get()
        self.UTR3p = magic.UTR3p.get()
        self.transcript_name = magic.transcript_name.get()
        self.construction_name = magic.construction_name.get()
        self.exon_number = magic.exon_number.get()
        self.exs_list = magic.exs_list

        self.values = {"output directory": self.outdir, "input file": self.seq, "constitutive exon vector 5 prime sequence": self.UTR5p, "constitutive exon vector 3 prime sequence": self.UTR3p, "gene name": self.gene_name, "transcript name": self.transcript_name, "construction name": self.construction_name, "exon number and sequence": self.exs_list}
        self.values_not_completed = [val_name for val_name, value in self.values.items() if not value]

        # generate files if all fields are completed
        if self.values_not_completed:
            message = "These fields are not defined: \n\n\t" + "\n\t".join(self.values_not_completed) + "\n\nDid you click 'Validate'?\n\nEXITING..."
            messagebox.showerror("ERROR. Files not generated", message)
            print("These fields are not defined: " + ", ".join(self.values_not_completed))
            return
            
        self.create_fasta_file()
        self.create_gtf_file()
    

    def create_fasta_file(self):
        self.seq_search_UTR5p_start = self.seq.lower().find(self.UTR5p.lower())
        self.seq_search_UTR3p_end = self.seq.lower().find(self.UTR3p.lower())

        self.insert_start = self.seq[self.seq_search_UTR5p_start+len(self.UTR5p):]
        self.insert_end = self.seq[:self.seq_search_UTR3p_end]

        self.seq_split = self.UTR5p + self.insert_start + self.insert_end + self.UTR3p
        
        with open(self.outdir+f"/{self.construction_name}.fasta", "w") as filout:
            filout.write(f'#MAGIC creation of artficial genome fasta file for {self.construction_name} construction used for minigene splicing assays and massively parallel sequencing\n')
            filout.write(f">{self.construction_name}\n")
            filout.write(self.seq_split)


    def create_gtf_file(self):
        seq_UTR5_start = self.seq_split.lower().find(self.UTR5p.lower()); seq_UTR5_end = seq_UTR5_start + len(self.UTR5p)
        seq_UTR3_start = self.seq_split.lower().find(self.UTR3p.lower()); seq_UTR3_end = seq_UTR3_start + len(self.UTR3p)
        strand = "+"
        
        with open (self.outdir+f"/{self.construction_name}.gtf", "w") as filout:
            filout.write(f'#MAGIC creation of artficial genome gtf file for {self.construction_name} construction used for minigene splicing assays and massively parallel sequencing\n')
            filout.write(f'{self.construction_name}\tMAGIC\ttranscript\t1\t{len(self.seq_split)}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}";\n')
            filout.write(f'{self.construction_name}\tMAGIC\tUTR\t{str(seq_UTR5_start+1)}\t{seq_UTR5_end}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}";\n')
            for ex in self.exs_list:
                ex_nb = ex[0]; ex_seq = ex[1]
                self.seq_search_start = self.seq_split.lower().find(ex_seq.lower())
                self.seq_search_end = self.seq_search_start + len(ex_seq)
                filout.write(f'{self.construction_name}\tMAGIC\tCDS\t{str(self.seq_search_start+1)}\t{str(self.seq_search_end)}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}"; exon_number "{ex_nb}";\n')
                filout.write(f'{self.construction_name}\tMAGIC\texon\t{str(self.seq_search_start+1)}\t{str(self.seq_search_end)}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}"; exon_number "{ex_nb}";\n')
            filout.write(f'{self.construction_name}\tMAGIC\tUTR\t{str(seq_UTR3_start+1)}\t{seq_UTR3_end}\t.\t{strand}\t.\tgene_name "{self.gene_name}"; transcript_id "{self.transcript_name}_{self.construction_name}";\n')



if __name__ == "__main__":
    magic = MAGIC_GUI()
    magic.mainloop()
    generate_files()