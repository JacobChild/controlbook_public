import VTOLParam as P
import hw16bodePlantCtrl as P16 
from control import bode, margin, mag2db
import matplotlib.pyplot as plt

# flag to define if using dB or absolute scale for M(omega)
dB_flag = P16.dB_flag

#Assign the transfer functions from the previous HWs
FtoH = P16.FtoH
TautoTheta = P16.TautoTheta
ThetatoZ = P16.ThetatoZ
C_pid_FtoH = P16.C_pid_FtoH
C_pid_TautoTheta = P16.C_pid_TautoTheta
C_pid_ThetatoZ = P16.C_pid_ThetatoZ

def MyBodeFunc(Plant,C_pid,title):
    
    if __name__ == '__main__':
        #calculate the phase and gain margins
        if dB_flag:
            gm, pm, Wcg, Wcp = margin(Plant*C_pid)
            gm = mag2db(gm)
            print("Inner Loop: ", "gm: ", gm, ", pm: ", pm, ", Wcg: ", Wcg, ", Wcp: ", Wcp)
        else:
            gm, pm, Wcg, Wcp = margin(Plant*C_pid)
            print("Inner Loop: ", "gm: ", gm, ", pm: ", pm, ", Wcg: ", Wcg, ", Wcp: ", Wcp)
            
        #display bode plots of transfer functions
        fig1 = plt.figure()
        
        #this makes two bode plots for open and closed loop
        bode(Plant * C_pid, dB=dB_flag,
            label='C(s)P(s) - Open-loop')
        bode(Plant*C_pid/(1+Plant*C_pid), dB=dB_flag,
            label=r'$\frac{P(s)C(s)}{1+P(s)C(s)}$ - Closed-loop')
        # now we can add lines to show where we calculated the GM and PM
        gm_line = fig1.axes[0].plot([Wcg, Wcg],
                                    plt.ylim(), 'k--', label='GM')
        gm_line[0].set_label('GM')
        fig1.axes[0].legend()
        fig1.axes[1].plot([Wcp, Wcp], plt.ylim(), 'b--', label='PM')
        plt.legend()

        # setting axis title
        fig1.axes[0].set_title(title +
                            'GM:'+str(round(gm, 2))+
                            ', PM:'+str(round(pm, 2)))

        plt.show()
        
        
MyBodeFunc(FtoH,C_pid_FtoH,"FtoH")
MyBodeFunc(TautoTheta,C_pid_TautoTheta,"TautoTheta")
MyBodeFunc(ThetatoZ,C_pid_ThetatoZ,"ThetatoZ")