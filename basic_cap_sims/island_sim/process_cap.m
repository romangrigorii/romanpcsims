dat = importdata('C:\Users\atrox\Desktop\Shared_Folder\sim_files_cap\capacitor_vals_noi.txt')

d = [.03 .05 .07 .09 .11 .13 .15];
c = dat/(1e-12)/5;

plot(d,c)
