clear; close all;

current_folder = pwd;
system(sprintf('cd %s', current_folder));
cmd = sprintf('CatGT -dir=%s -run=Exp -g=0 -t=0 -ob -ap -prb=0 -obx=0 -prb_fld -dest=%s -zerofillmax=0 -gblcar', current_folder, current_folder);
system(cmd);

% pause(1500);

movefile(fullfile(current_folder, 'extractObxData.m'), fullfile(current_folder, 'catgt_Exp_g0'));
cd('catgt_Exp_g0\');
extractObxData;

cd('../');
main_kilosort;

mainAutoCuration;
