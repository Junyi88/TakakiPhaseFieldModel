clear

clear;

Prefix0=['../exec/'];


for nt=0
    filename=[Prefix0 'Test_Phi_' num2str(nt) '.csv'];
    Phi=csvread(filename);

    filename=[Prefix0 'Test_Theta_' num2str(nt)  '.csv'];
    Theta=csvread(filename);

    figure(1);
    clf;
    subplot(1,2,1);
    surf(Phi,'EdgeColor','none');
    title(['\phi at nTime =' num2str(nt)]);
    view([0 90]);
    xlim([0 211]);
    ylim([0 211]);
    pbaspect([1 1 1]);
    

    subplot(1,2,2);
    surf(Theta,'EdgeColor','none');
    title(['\theta at nTime =' num2str(nt)]);
    view([0 90]);
    xlim([0 211]);
    ylim([0 211]);
    pbaspect([1 1 1]);
    
    pause(0.01);
end