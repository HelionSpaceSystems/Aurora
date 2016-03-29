# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 20:01:53 2016

@author: aliba
"""

import tkinter as tk
from tkinter.ttk import*

class Application(Frame):
    def __init__(self, master):
        Frame.__init__(self, master)
        self.grid()
        self.create_control_frame()


    def create_control_frame(self):
        
        def get_user_input():
            for results in self.inputs:
                print(results.get())
        
        self.inputs = []
        launch_labels = ['Launch Site: ', 'Launch Date: ', 'Payload Mass [kg]: ']
        atm_labels = ['Atmospheric Model: ', 'Wind Model: ']
        stage_labels = ['Vertical Ascent [s]: ', 'Coast Time [s]: ', 'Pitch Duration [s]: ', 'Burnout [s]: ', 'Engine Cutoff to Separation [s]: ']
        stage_tab_names = ['First Stage','Second Stage','Third Stage']
        tab_label_frames = ['Time Events','Aerodynamic Properties', 'Vehicle Properties']
        control_label_frame_names = ['Launch Properties', 'Desired Orbit','Atmospheric Models', 'Stage Properties']
        aero_labels = ['Lift Coeficient: ','Drag Coefficient: ','Referance Area [m^2]: ']
        vehicle_labels = ['Fuel Type: ', 'Oxidizer: ','Mixture Ratio: ',
                          'Payload Ratio: ','Exit Velocity [km/s]: ',
                          'Ratio of Specific Heats: ','Chamber Temperature [k]: ',
                          'Chamber Pressure [MPa]: ','Contraction Ratio: ',
                          'Design Altitude [km]: ','Structural Ratio: ','Max G: ',
                          'No. of Engines: ']
                          
        desiredOrbit_labels = ['Semimajor Axis [km]: ','Inclination [deg]: ',
                               'RAAN [deg]: ','Argument of Perigee [deg]: ',
                               'Angular Momentum [km^3/s^2]: ','Eccentricity: ',
                               'True Anomaly: ']
      
        for j in range(0,len(control_label_frame_names)):
            
            
            self.labelframe = LabelFrame(self, text = control_label_frame_names[j])
            self.labelframe.grid(row = j, padx = 5, pady = 5, sticky = 'NEWS')
            
            
            if control_label_frame_names[j] == 'Launch Properties':
                
                for k in range(0,len(launch_labels)):
                     self.label = Label(self.labelframe, text = launch_labels[k])
                     self.label.grid(row = k, padx = 5, pady = 5, sticky = 'NEWS')
                    
            if control_label_frame_names[j] == 'Desired Orbit':
                
                for k in range(0,len(desiredOrbit_labels)):
                     self.label = Label(self.labelframe, text = desiredOrbit_labels[k])
                     self.label.grid(row = k, padx = 5, pady = 5, sticky = 'NEWS')
                     
            if control_label_frame_names[j] == 'Atmospheric Models':
                
                for k in range(0,len(atm_labels)):
                     self.label = Label(self.labelframe, text = atm_labels[k])
                     self.label.grid(row = k, padx = 5, pady = 5, sticky = 'NEWS')
            
            if control_label_frame_names[j] == 'Stage Properties':
            
                self.notebook = Notebook(self.labelframe)
                for k in range(0, len(stage_tab_names)):
                    self.tab = Frame(self.notebook)
                    self.notebook.add(self.tab, text = stage_tab_names[k])
                    self.notebook.grid(row = k, column = 0, padx = 5, pady = 5, sticky = 'NEWS')
                    
                    for k in range(0, len(tab_label_frames)):
                        self.labelframe = LabelFrame(self.tab, text = tab_label_frames[k])
                        self.labelframe.grid(row = k, column = 0, padx = 5, pady = 5, sticky = 'NEWS')
                        
                        if tab_label_frames[k] == 'Time Events':
                            for k in range(0,len(stage_labels)):
                                self.label = Label(self.labelframe, text = stage_labels[k])
                                self.label.grid(row = k, column = 0, padx = 5, pady = 5, sticky = 'NEWS')
                                self.entry = Entry(self.labelframe)
                                self.entry.grid(row = k, column = 1, padx = 5, pady = 5, sticky = 'NEWS')
                                self.inputs.append(self.entry.get())
                        
                        elif tab_label_frames[k] == 'Aerodynamic Properties':
                            for k in range(0, len(aero_labels)):
                                    self.label = Label(self.labelframe, text = aero_labels[k])
                                    self.label.grid(row = k, column = 0, padx = 5, pady = 5, sticky = 'NEWS')
                                    self.entry = Entry(self.labelframe)
                                    self.entry.grid(row = k, column = 1, padx = 5, pady = 5, sticky = 'NEWS')
                                    self.inputs.append(self.entry)

                        elif tab_label_frames[k] == 'Vehicle Properties':
                            for k in range(0, len(vehicle_labels)):
                                    self.label = Label(self.labelframe, text = vehicle_labels[k])
                                    self.label.grid(row = k, column = 0, padx = 5, pady = 5, sticky = 'NEWS')
                                    self.entry = Entry(self.labelframe)
                                    self.entry.grid(row = k, column = 1, padx = 5, pady = 5, sticky = 'NEWS')
                                    #self.inputs.append(self.entry.get())
                                    #self.button = Button(self.labelframe,command = get_user_input).grid(row=k)
                                    
                                    

                            
                        
                       
                       
            
        
        
if __name__ == '__main__':
    root = tk.Tk()
    root.title("Aurora")
    root.geometry("1400x800")
    app = Application(root)
    root.mainloop()

