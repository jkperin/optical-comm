function m_file_list = get_m_file_list(base_dir)
%% Get list of all m files in current folder "base_dir" and all its subfolders
disp(base_dir)

files = dir([base_dir '/*']);
m_file_list = {};
for k = 1:length(files)
    if files(k).isdir
        if string(files(k).name) == "." || string(files(k).name) == ".." ||  contains(files(k).name, '.git')
            continue
        end
        m_file_list_sub = get_m_file_list([base_dir '\' files(k).name]);
        m_file_list = [m_file_list m_file_list_sub];
    elseif contains(files(k).name, '.m')
        full_path = [files(k).folder '\' files(k).name];
        disp(full_path)
        m_file_list = [m_file_list full_path];
    else
        continue
    end
end