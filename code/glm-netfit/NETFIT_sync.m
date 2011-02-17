%SCRIPT SYNCHRONIZATION POINT FOR PARALLEL EXECS

if(N_proc==1) return; end

%INIT SYNC SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist('sync_step')~=1) sync_step=1; end
fsname=[sync_name,sprintf('-lock%i.mat',id_proc)];
fcname=[sync_name,sprintf('-release%i.mat',id_proc)];
fbname=[sync_name,sprintf('-buffer%i.mat',id_proc)];
if(exist('syncdata')~=1) syncdata={}; end
delete(fcname);
save(fsname,'N_proc','id_proc','sync_step');
pause(5)
delete(fbname);
save(fbname,syncdata{:});
pause(5);

%SYNCHRONIZE ARRIVALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flg=1;
while(flg)
  flg=0;
  for k=1:N_proc
    fsname_tmp=[sync_name,sprintf('-lock%i.mat',k)];      
    if(exist(fsname_tmp)~=2) flg=1; end
  end

  fprintf('Sync lock waiting...\n');
  if(flg) pause(15); end            %wait for others, then check again/exit
end
save(fcname,'N_proc','id_proc','sync_step');
pause(5)

flg=1;
while(flg)
  flg=0;
  for k=1:N_proc
    fcname_tmp=[sync_name,sprintf('-release%i.mat',k)]; 
    if(exist(fcname_tmp)~=2) flg=1; end
  end

  fprintf('Sync release waiting...\n');
  if(flg) pause(15); end            %wait for others, then check again/exit
end
delete(fsname);


%SYNCHRONIZE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for isync=1:length(syncdata)
    var=syncdata{isync};
    tmp=eval(var);
    for ksync=1:N_proc
        fbname=[sync_name,sprintf('-buffer%i.mat',ksync)];
        load(fbname,var);
        for lsync=nrange{ksync}
            if(lsync==0) continue; end
            tmp{lsync}=eval([var,sprintf('{%i}',lsync)]);
        end
    end
    eval([var,'=tmp;']);
end

fprintf('Barrier========================\n');
pause(120);

clear tmp flg fsname fcname fbname k fsname_tmp fcname_tmp var isync ksync lsync