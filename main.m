for isCollaborate = 0:2
    main_distribute;
    if isCollaborate == 0
        save('../autonomous.mat');
    elseif isCollaborate == 2
        save('../collaborate_feedin.mat');
    else
        save('../collaborate.mat');
    end
end
main_handle_171013_v2
