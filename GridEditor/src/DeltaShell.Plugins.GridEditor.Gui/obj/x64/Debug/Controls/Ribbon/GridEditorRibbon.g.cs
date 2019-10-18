﻿#pragma checksum "..\..\..\..\..\Controls\Ribbon\GridEditorRibbon.xaml" "{8829d00f-11b8-4213-878b-770e8597ac16}" "9F0858AC4BD046344C5867EFC68833FDF8C65B80DA0C9E73C5B8567C4C474FA0"
//------------------------------------------------------------------------------
// <auto-generated>
//     This code was generated by a tool.
//     Runtime Version:4.0.30319.42000
//
//     Changes to this file may cause incorrect behavior and will be lost if
//     the code is regenerated.
// </auto-generated>
//------------------------------------------------------------------------------

using DelftTools.Controls.Wpf.Extensions;
using DelftTools.Controls.Wpf.ValueConverters;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon;
using DeltaShell.Plugins.GridEditor.Helpers;
using Fluent;
using NetTopologySuite.Extensions.Grids;
using System;
using System.Diagnostics;
using System.Windows;
using System.Windows.Automation;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Ink;
using System.Windows.Input;
using System.Windows.Markup;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Media.Effects;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Windows.Media.TextFormatting;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Shell;


namespace DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon {
    
    
    internal partial class GridEditorRibbon : System.Windows.Controls.UserControl, System.Windows.Markup.IComponentConnector {
        
        
        #line 42 "..\..\..\..\..\Controls\Ribbon\GridEditorRibbon.xaml"
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1823:AvoidUnusedPrivateFields")]
        private Fluent.Ribbon RibbonControl;
        
        #line default
        #line hidden
        
        
        #line 44 "..\..\..\..\..\Controls\Ribbon\GridEditorRibbon.xaml"
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1823:AvoidUnusedPrivateFields")]
        internal Fluent.RibbonContextualTabGroup geospatialContextualGroup;
        
        #line default
        #line hidden
        
        
        #line 47 "..\..\..\..\..\Controls\Ribbon\GridEditorRibbon.xaml"
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1823:AvoidUnusedPrivateFields")]
        internal Fluent.RibbonTabItem GridEditorTab;
        
        #line default
        #line hidden
        
        
        #line 52 "..\..\..\..\..\Controls\Ribbon\GridEditorRibbon.xaml"
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1823:AvoidUnusedPrivateFields")]
        internal DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.GridEditorRibbonViewModel GridEditingGroupViewModel;
        
        #line default
        #line hidden
        
        private bool _contentLoaded;
        
        /// <summary>
        /// InitializeComponent
        /// </summary>
        [System.Diagnostics.DebuggerNonUserCodeAttribute()]
        [System.CodeDom.Compiler.GeneratedCodeAttribute("PresentationBuildTasks", "4.0.0.0")]
        public void InitializeComponent() {
            if (_contentLoaded) {
                return;
            }
            _contentLoaded = true;
            System.Uri resourceLocater = new System.Uri("/DeltaShell.Plugins.GridEditor.Gui;component/controls/ribbon/grideditorribbon.xam" +
                    "l", System.UriKind.Relative);
            
            #line 1 "..\..\..\..\..\Controls\Ribbon\GridEditorRibbon.xaml"
            System.Windows.Application.LoadComponent(this, resourceLocater);
            
            #line default
            #line hidden
        }
        
        [System.Diagnostics.DebuggerNonUserCodeAttribute()]
        [System.CodeDom.Compiler.GeneratedCodeAttribute("PresentationBuildTasks", "4.0.0.0")]
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1811:AvoidUncalledPrivateCode")]
        internal System.Delegate _CreateDelegate(System.Type delegateType, string handler) {
            return System.Delegate.CreateDelegate(delegateType, this, handler);
        }
        
        [System.Diagnostics.DebuggerNonUserCodeAttribute()]
        [System.CodeDom.Compiler.GeneratedCodeAttribute("PresentationBuildTasks", "4.0.0.0")]
        [System.ComponentModel.EditorBrowsableAttribute(System.ComponentModel.EditorBrowsableState.Never)]
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Design", "CA1033:InterfaceMethodsShouldBeCallableByChildTypes")]
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Maintainability", "CA1502:AvoidExcessiveComplexity")]
        [System.Diagnostics.CodeAnalysis.SuppressMessageAttribute("Microsoft.Performance", "CA1800:DoNotCastUnnecessarily")]
        void System.Windows.Markup.IComponentConnector.Connect(int connectionId, object target) {
            switch (connectionId)
            {
            case 1:
            this.RibbonControl = ((Fluent.Ribbon)(target));
            return;
            case 2:
            this.geospatialContextualGroup = ((Fluent.RibbonContextualTabGroup)(target));
            return;
            case 3:
            this.GridEditorTab = ((Fluent.RibbonTabItem)(target));
            return;
            case 4:
            this.GridEditingGroupViewModel = ((DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.GridEditorRibbonViewModel)(target));
            return;
            }
            this._contentLoaded = true;
        }
    }
}

