


clf
subplot(4,4,[1 5 9])
imagesc(tot_ndat.active2nrem(kidx_active2nrem==4,:));
title('active2nrem')
caxis([-0. 6])
subplot(4,4,13)
shadedErrorBar(1:1201,mean(tot_ndat.active2nrem(kidx_active2nrem==4,:)),std(tot_ndat.active2nrem(kidx_active2nrem==4,:))/sqrt(size(tot_ndat.active2nrem(kidx_active2nrem==4,:),1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[2 6 10])
imagesc(tot_ndat.nrem2active(kidx_nrem2active==3,:));
title('nrem2active')
caxis([0.5 6])
subplot(4,4,14)
shadedErrorBar(1:1201,mean(tot_ndat.nrem2active(kidx_nrem2active==3,:)),std(tot_ndat.nrem2active(kidx_nrem2active==3,:))/sqrt(size(tot_ndat.nrem2active(kidx_nrem2active==3,:),1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[3 7 11])
imagesc(tot_ndat.nrem2rem(kidx_nrem2rem==3,:))
title('nrem2rem')
caxis([0 6])
subplot(4,4,15)
shadedErrorBar(1:1201,mean(tot_ndat.nrem2rem(kidx_nrem2rem==3,:)),std(tot_ndat.nrem2rem(kidx_nrem2rem==3,:))/sqrt(size(tot_ndat.nrem2rem(kidx_nrem2rem==3,:),1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[4 8 12])
imagesc(tot_ndat.rem2nrem(kidx_rem2nrem==2,:));
caxis([0 6])
title('rem2nrem')
subplot(4,4,16)
shadedErrorBar(1:1201,mean(tot_ndat.rem2nrem(kidx_rem2nrem==2,:)),std(tot_ndat.rem2nrem(kidx_rem2nrem==2,:))/sqrt(size(tot_ndat.rem2nrem(kidx_rem2nrem==2,:),1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])




clf
subplot(4,4,[1 5 9])
imagesc(tot1215.active2nrem(kidx_active2nrem==2,:));
title('active2nrem')
caxis([-0.5 6])
subplot(4,4,13)
shadedErrorBar(1:1201,mean(tot1215.active2nrem(kidx_active2nrem==2,:)),std(tot1215.active2nrem(kidx_active2nrem==2,:))/sqrt(size(tot_ndat.active2nrem(kidx_active2nrem==4,:),1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[2 6 10])
imagesc(tot_ndat.nrem2active(kidx_nrem2active==3,:));
title('nrem2active')
caxis([0.5 6])
subplot(4,4,14)
shadedErrorBar(1:1201,mean(tot_ndat.nrem2active(kidx_nrem2active==3,:)),std(tot_ndat.nrem2active(kidx_nrem2active==3,:))/sqrt(size(tot_ndat.nrem2active(kidx_nrem2active==3,:),1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[3 7 11])
imagesc(tot_ndat.nrem2rem(kidx_nrem2rem==1,:))
title('nrem2rem')
caxis([0 6])
subplot(4,4,15)
shadedErrorBar(1:1201,mean(tot_ndat.nrem2rem(kidx_nrem2rem==1,:)),std(tot_ndat.nrem2rem(kidx_nrem2rem==1,:))/sqrt(size(tot_ndat.nrem2rem(kidx_nrem2rem==1,:),1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[4 8 12])
imagesc(tot_ndat.rem2nrem);
caxis([0 6])
title('rem2nrem')
subplot(4,4,16)
shadedErrorBar(1:1201,mean(tot_ndat.rem2nrem),std(tot_ndat.rem2nrem)/sqrt(size(tot_ndat.rem2nrem,1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])











clf
subplot(4,4,[1 5 9])
imagesc(tot_ndat.active2nrem(kidx_active2nrem==4,:));
title('active2nrem')
caxis([-0. 6])
subplot(4,4,13)
shadedErrorBar(1:1201,mean(tot_ndat.active2nrem(kidx_active2nrem==4,:)),std(tot_ndat.active2nrem(kidx_active2nrem==4,:))/sqrt(size(tot_ndat.active2nrem(kidx_active2nrem==4,:),1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])


clf
subplot(4,4,[2 6 10])
tot_ndat.rem2qwake = unique(tot_ndat.rem2qwake,'rows');
imagesc(tot_ndat.rem2qwake);
title('rem 2 qwake')
caxis([0 5])
subplot(4,4,14)
shadedErrorBar(1:1201,mean(tot_ndat.rem2qwake),std(tot_ndat.rem2qwake)/sqrt(size(tot_ndat.rem2qwake,1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[3 7 11])

imagesc(tot_ndat.qwake2awake(kidx_qwake2awake==3,:))
title('quiet 2 awake')
caxis([0 5])
subplot(4,4,15)
shadedErrorBar(1:1201,mean(tot_ndat.qwake2awake(kidx_qwake2awake==3,:)),std(tot_ndat.qwake2awake(kidx_qwake2awake==3,:))/sqrt(size(tot_ndat.qwake2awake(kidx_qwake2awake==3,:),1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[4 8 12])
imagesc(tot_ndat.awake2qwake(kidx_awake2qwake==2,:));
caxis([0 5])
title('active 2 awake')
subplot(4,4,16)
shadedErrorBar(1:1201,mean(tot_ndat.awake2qwake(kidx_awake2qwake==2,:)),std(tot_ndat.awake2qwake(kidx_awake2qwake==2,:))/sqrt(size(tot_ndat.awake2qwake(kidx_awake2qwake==2,:),1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])


exportgraphics(gcf,fullfile(dirpath,['zt0_3_active_quiet_transitions.pdf']),'Resolution',300,'ContentType','vector')



clf

subplot(4,4,[1 5 9])
tot_ndat.qwake2rem = unique(tot_ndat.qwake2rem,'rows');
imagesc(tot_ndat.qwake2rem);
title('qwake 2 rem')
caxis([-0. 6])
subplot(4,4,13)
plot(tot_ndat.qwake2rem)
%shadedErrorBar(1:1201,mean(tot_ndat.qwake2rem),std(tot_ndat.qwake2rem)/sqrt(size(tot_ndat.qwake2rem,1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])



subplot(4,4,[2 6 10])
tot_ndat.rem2qwake = unique(tot_ndat.rem2qwake,'rows');
imagesc(tot_ndat.rem2qwake);
title('rem 2 qwake')
caxis([0 5])
subplot(4,4,14)
shadedErrorBar(1:1201,mean(tot_ndat.rem2qwake),std(tot_ndat.rem2qwake)/sqrt(size(tot_ndat.rem2qwake,1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[3 7 11])
%tot_ndat.qwake2awake = unique(tot_ndat.qwake2awake,'rows');
imagesc(tot_ndat.qwake2awake(kidx_qwake2awake==1,:))
title('quiet 2 awake')
caxis([0 5])
subplot(4,4,15)
shadedErrorBar(1:1201,mean(tot_ndat.qwake2awake(kidx_qwake2awake==1,:)),std(tot_ndat.qwake2awake(kidx_qwake2awake==1,:))/sqrt(size(tot_ndat.qwake2awake(kidx_qwake2awake==1,:),1)));
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])

subplot(4,4,[4 8 12])
tot_ndat.awake2qwake = unique(tot_ndat.awake2qwake,'rows');
imagesc(tot_ndat.awake2qwake(kidx_awake2qwake==2,:));
caxis([0 5])
title('active 2 qwake')
subplot(4,4,16)
shadedErrorBar(1:1201,mean(tot_ndat.awake2qwake(kidx_awake2qwake==2,:)),std(tot_ndat.awake2qwake(kidx_awake2qwake==2,:))/sqrt(size(tot_ndat.awake2qwake(kidx_awake2qwake==2,:),1)))
xticks([0 600 1200])
xticklabels({'-24' '0' '24'})
xlabel('Time (s)')
xlim([0 1200])





colorbar



clear tot_ndat
for k=1:numel(fn)
    subplot(1,4,k)
    tot_ndat.(fn{k}) = tot1215.(fn{k}) - tot1215.(fn{k})(:,1);
    plot(mean(tot1215.(fn{k})),'LineWidth',1)
    title(fn{k})
    xlim([0 1200])
    xticks([0 600 1200])
    grid on
    ax = gca;
    ax.GridColor = [1, 0, 0];
    ax.GridAlpha = 0.75;
    ax.YGrid = 'off';
    hold all
end







