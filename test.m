System = 'data/S0_50N0_10SIM';

% Load 
load([System '/fov1_track.mat']);
% Calculate MSD
knownMSD = calcMSD(Tracks,1,0);
load([System '/movie4QMData.mat'])
QMMSD = correctedMSDs;
preTracks = Tracks;
preTracks(:,3:4) = [];
preMSD = calcMSD(preTracks,1,0);
%Plot MSDs

fig1 = figure();
whitebg(fig1,[1,1,1])
scatter(1:size(knownMSD,1),knownMSD(:,1),'k')
hold on
errorbar(1:size(preMSD,1),preMSD(:,1),preMSD(:,2)/sqrt(size(preMSD,1)),'o','Color','g')
hold on
errorbar(1:size(MSDs,1),MSDs(:,1),MSDs(:,2)/sqrt(size(MSDs,1)),'o','Color','b')
hold on
errorbar(1:size(QMMSD,1),QMMSD(:,1),QMMSD(:,2)/sqrt(size(QMMSD,1)),'o','Color','r')
ylim([min(QMMSD(2:end,1))*0.8,max(QMMSD(:,1))*1.5])
box on
set(gca,'xscale','log','yscale','log')
xlabel('\tau [frame]');
ylabel('\langledR^{2}\rangle [nm^{2}]');
legend('known','preMSD','4QM','4QM-corrected','Location','northwest')
legend('boxoff')
