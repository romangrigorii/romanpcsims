dat = importdata('C:\Users\atrox\Desktop\Shared_Folder\sim_files_cap2\capacitor_vals_i.txt')

d = [.03 .05 .07 .09 .11 .13 .15];
c = dat/(1e-12)/5;

hold on
a = plot(d,c,'o');
a.Color = [1 0 0];
a = plot(d,c);
a.Color = [1 0 0];

axis([0 .16 0 .3])