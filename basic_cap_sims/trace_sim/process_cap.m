source =  'C:\Users\atrox\Desktop\Shared_Folder\basic_cap_sims\trace_sim\UNV\';
d = round([.02,.04,.06,.08,.1,.12,.14,.16,.18,.2,.22,.24,.26,.28,.3]*1000);
c = [];
for i = 1:length(d)
    if length(num2str(d(i))) == 2
        dat = importdata(strcat(source,'Mesh_1_00',num2str(d(i)),'\mesh_1_00',num2str(d(i)),'_cap'));
    else
        dat = importdata(strcat(source,'Mesh_1_0',num2str(d(i)),'\mesh_1_0',num2str(d(i)),'_cap'));
    end
    c(i) = sum(dat(2:end)/(1e-12)/5);
end

