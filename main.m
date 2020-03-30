for isCollaborate = 0:2
    main_distribute;
    if isCollaborate == 0
        save('data/autonomous.mat');
    elseif isCollaborate == 2
        save('data/collaborate_feedin.mat');
    else
        save('data/collaborate.mat');
    end
end
main_handle
