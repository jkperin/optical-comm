%% Check dependencies

base_dir = pwd;

%%
m_file_list = get_m_file_list([base_dir '/apd']);
dep = dependencies.toolboxDependencyAnalysis(m_file_list)

%%
m_file_list = get_m_file_list([base_dir '/coherent']);
dep = dependencies.toolboxDependencyAnalysis(m_file_list)

%%
m_file_list = get_m_file_list([base_dir '/edfa']);
dep = dependencies.toolboxDependencyAnalysis(m_file_list)

%%
m_file_list = get_m_file_list([base_dir '/f']);
dep = dependencies.toolboxDependencyAnalysis(m_file_list)

%%
m_file_list = get_m_file_list([base_dir '/mpam']);
dep = dependencies.toolboxDependencyAnalysis(m_file_list)

%%
m_file_list = get_m_file_list([base_dir '/ofdm']);
dep = dependencies.toolboxDependencyAnalysis(m_file_list)

%%
m_file_list = get_m_file_list([base_dir '/soa']);
dep = dependencies.toolboxDependencyAnalysis(m_file_list)



files = {'C:\Users\Joe\Dropbox\research\codes\coherent\f\ber_coherent_awgn.m',...
    'C:\Users\Joe\Dropbox\research\codes\coherent\ber_coherent_dsp.m',...
    'C:\Users\Joe\Dropbox\research\codes\coherent\DSP\fse_cma.m',...
    'C:\Users\Joe\Dropbox\research\codes\coherent\QPSK_DSP_BER_qsub.m'}
    