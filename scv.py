import flet as ft
from flet.matplotlib_chart import MatplotlibChart
import numpy as np
import pandas as pd
import matplotlib as mlt
import matplotlib.pyplot as plt
#import matplotlib.backends.backend_svg
#import matplotlib.backends.backend_tkagg
import seaborn as sns
import pickle

#-------------------------------------------------
l_xy=['PCA', 'UMAP', 'Harmony', 'tSNE']  #scatter plot types, for choosing plt method 
l_violin=['Sample', 'Cell']  #dropdown options if violin plot is selected

#-------------------------------------------------
def main(page):
    #misc
    page.title='scRNASeq Viewer'    

    #file pick event
    def evt_file(e):
        #load data
        with open(e.files[0].path, 'rb') as f: evt_file.obj=pickle.load(f)
        r0_txt.value=e.files[0].path
        #dropdown: plot
        l_plot=list(evt_file.obj['plot'])
        r2_plot.options=[ft.dropdown.Option(i) for i in l_plot]
        #dropdown: feature
        l_feature=list(evt_file.obj['feature'])
        r2_feature.options=[ft.dropdown.Option(i) for i in l_feature]
        #update
        page.update()
        return
    
    #plot type select event
    def evt_plot(e):
        #reset
        r2_feature.value=''
        r2_feature.disabled=False
        r2_value.options=[]
        r2_value.hint_text=''
        r2_value.disabled=True 
        r2_gene.value=''
        r2_gene.disabled=True
        #change 
        l_feature=list(evt_file.obj['feature'])
        if r2_plot.value=='Violin': l_feature=[i for i in l_feature if i in l_violin]
        r2_feature.options=[ft.dropdown.Option(i) for i in l_feature]
        #update
        page.update()
        return

    #feature select event
    def evt_feature(e):
        #reset
        r2_value.options=[] 
        r2_value.hint_text=''
        r2_value.disabled=True 
        r2_gene.value=''
        r2_gene.disabled=True
        #change
        if r2_feature.value=='Sample':
            l_sample=evt_file.obj['feature']['Sample']['value']
            l_sample=['All Samples']+l_sample    
            r2_value.options=[ft.dropdown.Option(i) for i in l_sample]
            r2_value.hint_text='Choose sample'
            r2_value.disabled=False
        if r2_feature.value=='Cell':
            l_cell=evt_file.obj['feature']['Cell']['value']
            if r2_plot.value!='Violin': l_cell=['All Cells']+l_cell
            r2_value.options=[ft.dropdown.Option(i) for i in l_cell]
            r2_value.hint_text='Choose cell'
            r2_value.disabled=False
        if r2_feature.value=='Gene':
            r2_gene.disabled=False
        #update
        if r2_plot.value=='Violin': r2_gene.disabled=False
        page.update()
        return
        
    #submit event
    def plt_blank():
        fig, ax=plt.subplots()
        plt.axis('off')
        return fig

    def plt_violin(dfi, title=None, cmap=None, key=None):
        #plot
        fig, ax=plt.subplots()
        sns.despine()
        ax=sns.violinplot(x='col', y='gene', data=dfi, palette=cmap, density_norm='width', linewidth=0.5, hue='col', legend=False)
        #adjust
        plt.title(title, weight='semibold', pad=10)
        plt.xlabel('')
        plt.ylabel('Normalized Counts', weight='semibold')
        if key=='sample': plt.xticks(rotation=30, ha='right', weight='semibold')
        if key=='cell': plt.xticks(rotation=0, weight='semibold')
        plt.tight_layout()
        return fig

    def plt_all(dfi, title=None, s=5, cmap='tab20'):
        #plot
        fig, ax=plt.subplots()
        plt.axis('off')
        ax=sns.scatterplot(x='x', y='y', data=dfi, hue='col', palette=cmap, s=s, alpha=0.7)
        plt.title(title, weight='semibold', pad=10)
        plt.legend(loc=(1.01, 0.1), frameon=False, markerscale=1, prop={'weight': 'semibold'})
        plt.tight_layout()
        return fig

    def plt_single(dfi, title=None, s=5, cmap=['#ff0000', '#888888']):
        #plot
        fig, ax=plt.subplots()
        plt.axis('off')
        ax=sns.scatterplot(x='x', y='y', data=dfi, hue='col', palette=cmap, s=s, alpha=0.7)
        plt.title(title, weight='semibold', pad=10)
        plt.legend(loc=(1.01, 0.1), frameon=False, markerscale=1, prop={'weight': 'semibold'})
        plt.tight_layout()
        return fig

    def plt_gene(dfi, title=None, cmap='BuPu', s=5):
        #prep
        vmax=dfi['col'].max()
        norm=plt.Normalize(-0.8, vmax)
        l_tick=np.arange(-0.8, vmax, 1.0).tolist()
        #plot
        fig, ax=plt.subplots()
        plt.axis('off')
        ax=sns.scatterplot(x='x', y='y', data=dfi, hue='col', palette=cmap, s=s, alpha=0.7, hue_norm=norm)
        plt.title(title, weight='semibold', pad=10)
        plt.legend(loc=(1.01, 0.1), frameon=False, markerscale=1, prop={'weight': 'semibold'})
        plt.tight_layout()
        return fig

    def evt_sub(e):
        check=True
        #param
        xy, feature, value, gene=r2_plot.value, r2_feature.value, r2_value.value, r2_gene.value
        gene=gene.strip()
        df=evt_file.obj['data']
        s=df.shape[0]/4000
        #cmap
        try: cmap_sample=evt_file.obj['feature']['Sample']['cmap']
        except Exception: cmap_sample='tab20' 
        try: cmap_cell=evt_file.obj['feature']['Cell']['cmap']
        except Exception: cmap_cell='tab20' 
        #violin - sample
        if (r2_plot.value=='Violin') and (r2_feature.value=='Sample'):
            #data
            df=df.reindex([gene, 'sample', 'cell'], axis=1).dropna(axis=1)
            if df.shape[1]==3: 
                if r2_value.value!='All Samples': df=df.query('sample==@r2_value.value')
                df=df.drop('sample', axis=1)
                df.columns=['gene', 'col']
                #plot
                title=f'{r2_gene.value} ({r2_value.value})'
                fig=plt_violin(df, title=title, cmap=cmap_cell, key='sample')
                check=False 
        #violin - cell, but no sample col
        if (r2_plot.value=='Violin') and (r2_feature.value=='Cell') and ('sample' not in df.columns.tolist()): 
            fig=plt_blank()
            check=False
        #violin - cell
        if (r2_plot.value=='Violin') and (r2_feature.value=='Cell') and ('sample' in df.columns.tolist()): 
            #data
            df=df.reindex([gene, 'sample', 'cell'], axis=1).dropna(axis=1)
            if df.shape[1]==3:
                df=df.query('cell==@r2_value.value')
                df=df.drop('cell', axis=1)
                df.columns=['gene', 'col']
                #plot
                title=f'{r2_gene.value} ({r2_value.value})'
                fig=plt_violin(df, title=title, cmap=cmap_sample, key='cell')
                check=False 
        #scatter - all sample
        if (r2_plot.value in l_xy) and (r2_feature.value=='Sample') and (r2_value.value=='All Samples'):
            #data
            df_xy=evt_file.obj['plot'][r2_plot.value]
            df=df.merge(df_xy, left_index=True, right_index=True).reindex(['x', 'y', 'sample'], axis=1)
            df=df.dropna(axis=1)
            if df.shape[1]==3:
                df.columns=['x', 'y', 'col']
                #plot
                title=f'All Samples ({r2_plot.value})'
                fig=plt_all(df, title=title, s=s, cmap=cmap_sample)
                check=False
        #scatter - all cell
        if (r2_plot.value in l_xy) and (r2_feature.value=='Cell') and (r2_value.value=='All Cells'):
            #data
            df_xy=evt_file.obj['plot'][r2_plot.value]
            df=df.merge(df_xy, left_index=True, right_index=True).reindex(['x', 'y', 'cell'], axis=1)
            df=df.dropna(axis=1)
            if df.shape[1]==3:
                df.columns=['x', 'y', 'col']
                #plot
                title=f'All Cells ({r2_plot.value})'
                fig=plt_all(df, title=title, s=s, cmap=cmap_cell)
                check=False
        #scatter - single sample
        if (r2_plot.value in l_xy) and (r2_feature.value=='Sample') and (r2_value.value!='All Samples'):
            #data
            df_xy=evt_file.obj['plot'][r2_plot.value]
            df=df.merge(df_xy, left_index=True, right_index=True).reindex(['x', 'y', 'sample'], axis=1)
            df=df.dropna(axis=1)
            if df.shape[1]==3:
                #prep
                df['tmp']='Others'
                df.loc[df['sample']==r2_value.value, ['tmp']]=r2_value.value
                df['tmp']=pd.Categorical(df['tmp'], categories=[r2_value.value, 'Others'], ordered=True)
                df.columns=['x', 'y', '_', 'col']
                #plot
                fig=plt_single(df, title=r2_value.value, s=s)
                check=False
        #scatter - single cell
        if (r2_plot.value in l_xy) and (r2_feature.value=='Cell') and (r2_value.value!='All Cells'):
            #data
            df_xy=evt_file.obj['plot'][r2_plot.value]
            df=df.merge(df_xy, left_index=True, right_index=True).reindex(['x', 'y', 'cell'], axis=1)
            df=df.dropna(axis=1)
            if df.shape[1]==3:
                #prep
                df['tmp']='Others'
                df.loc[df['cell']==r2_value.value, ['tmp']]=r2_value.value
                df['tmp']=pd.Categorical(df['tmp'], categories=[r2_value.value, 'Others'], ordered=True)
                df.columns=['x', 'y', '_', 'col']
                #plot
                fig=plt_single(df, title=r2_value.value, s=s)
                check=False
        #scatter - gene
        if (r2_plot.value in l_xy) and (r2_feature.value=='Gene'):
            #data
            df_xy=evt_file.obj['plot'][r2_plot.value]
            df=df.merge(df_xy, left_index=True, right_index=True).reindex(['x', 'y', gene], axis=1)
            df=df.dropna(axis=1)
            if df.shape[1]==3:
                df.columns=['x', 'y', 'col']
                #plot
                title=f'{gene} ({r2_plot.value})'
                fig=plt_gene(df, title=title, s=s)
                check=False
        #update
        if check: fig=plt_blank()
        r4_plot.figure=fig
        r4_plot.visible=True
        page.update()
        return

    #-------------------------------------------------
    #row_0: select file 
    r0_file=ft.FilePicker(on_result=evt_file)
    page.overlay.append(r0_file)
    r0_btn=ft.ElevatedButton("Choose a file...", on_click=lambda _: r0_file.pick_files(allow_multiple=False, allowed_extensions=['pkl']))
    r0_txt=ft.Text(value='No file selected')
    r0_container=ft.Row([r0_btn, r0_txt])

    #row_1: text
    r1_txt=ft.Text('Plot Settings:', size=15)

    #row_2: plot settings
    r2_plot=ft.Dropdown(options=[], hint_text='Plot Type', on_change=evt_plot) 
    r2_feature=ft.Dropdown(options=[], hint_text='Feature', on_change=evt_feature)
    r2_feature.disabled=True
    r2_value=ft.Dropdown(options=[], hint_text='') 
    r2_value.disabled=True
    r2_gene=ft.TextField(label='Gene Name')
    r2_gene.disabled=True
    r2_container=ft.Row([r2_plot, r2_feature, r2_value, r2_gene])

    #row_3: submit
    r3_sub=ft.ElevatedButton('Submit', on_click=evt_sub)

    #row_4: plot
    fig=plt_blank()
    r4_plot=MatplotlibChart(fig, expand=True)
    r4_plot.visible=False

    #update
    r0_div=ft.Divider(height=10, color='#66888888')
    r3_div=ft.Divider(height=10, color='#66888888')
    page.add(r0_container, r0_div, r1_txt, r2_container, r3_sub, r3_div, r4_plot)
    return

#-------------------------------------------------
ft.app(target=main)
