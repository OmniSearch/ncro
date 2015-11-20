<!-- Tufts VUE 3.2.2 concept-map (ncRO-UseCase-1-step-1.vue) 2015-06-08 -->
<!-- Tufts VUE: http://vue.tufts.edu/ -->
<!-- Do Not Remove: VUE mapping @version(1.1) jar:file:/C:/Program%20Files/VUE/VUE.jar!/tufts/vue/resources/lw_mapping_1_1.xml -->
<!-- Do Not Remove: Saved date Mon Jun 08 23:20:48 CDT 2015 by Jingshan Huang on platform Windows 7 6.1 in JVM 1.7.0_21-b11 -->
<!-- Do Not Remove: Saving version @(#)VUE: built May 23 2013 at 2146 by tomadm on Linux 2.6.18-348.2.1.el5 i386 JVM 1.7.0_21-b11(bits=32) -->
<?xml version="1.0" encoding="US-ASCII"?>
<LW-MAP xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:noNamespaceSchemaLocation="none" ID="0"
    label="ncRO-UseCase-1-step-1.vue" created="1431734769931" x="0.0"
    y="0.0" width="1.4E-45" height="1.4E-45" strokeWidth="0.0" autoSized="false">
    <resource referenceCreated="1433823648228" size="38397"
        spec="C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design\ncRO-UseCase-1-step-1.vue"
        type="1" xsi:type="URLResource">
        <title>ncRO-UseCase-1-step-1.vue</title>
        <property key="File" value="C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design\ncRO-UseCase-1-step-1.vue"/>
    </resource>
    <fillColor>#FFFFFF</fillColor>
    <strokeColor>#404040</strokeColor>
    <textColor>#000000</textColor>
    <font>SansSerif-plain-14</font>
    <URIString>http://vue.tufts.edu/rdf/resource/5a6f927f0a147a5701bb3ae4d29d2d28</URIString>
    <child ID="30" label="SO: miRNA" layerID="1" created="1431736299454"
        x="721.4205" y="324.5058" width="73.0" height="23.0"
        strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#9DDB53</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f929d0a147a5701bb3ae4fb5fe718</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="46" label="SO: gene" layerID="1" created="1431740614381"
        x="973.59924" y="326.1799" width="62.0" height="23.0"
        strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#9DDB53</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b00a147a5701bb3ae44cf956bd</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="47" label="ncro: miRNA_target_gene" layerID="1"
        created="1431740720414" x="620.79315" y="523.52997"
        width="153.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b00a147a5701bb3ae467b66d03</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="48" label="is_a" layerID="1" created="1431740743700"
        x="772.9625" y="348.67993" width="227.18457" height="189.98761"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b10a147a5701bb3ae465962ea4</URIString>
        <point1 x="773.4625" y="538.16754"/>
        <point2 x="999.6471" y="349.17993"/>
        <ID1 xsi:type="node">47</ID1>
        <ID2 xsi:type="node">46</ID2>
        <ctrlPoint0 x="915.74084" y="544.02826" xsi:type="point"/>
    </child>
    <child ID="59" label="ncro: miRNA_and_target_gene_binding"
        layerID="1" created="1433381847353" x="627.717" y="426.31097"
        width="229.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc3bc18dc0a80c610191f5df11a59903</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="60" label="ncro: is_about_miRNA_target_gene" layerID="1"
        created="1433381966836" x="645.01636" y="448.81104"
        width="171.0" height="75.218994" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc3bc18ec0a80c610191f5dfe25d1dd6</URIString>
        <point1 x="740.9772" y="449.31104"/>
        <point2 x="705.6995" y="523.53"/>
        <ID1 xsi:type="node">59</ID1>
        <ID2 xsi:type="node">47</ID2>
        <ctrlPoint0 x="737.6944" y="479.76086" xsi:type="point"/>
    </child>
    <child ID="61" label="ncro: is_about_miRNA" layerID="1"
        created="1433382183779" x="696.5687" y="347.00586" width="107.0"
        height="79.805176" strokeWidth="1.0" autoSized="false"
        controlCount="0" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc3cec2dc0a80c610191f5df3bd3b19f</URIString>
        <point1 x="743.99084" y="426.31104"/>
        <point2 x="756.1466" y="347.50586"/>
        <ID1 xsi:type="node">59</ID1>
        <ID2 xsi:type="node">30</ID2>
    </child>
    <child ID="63" label="ncro: human_miRNA" layerID="1"
        created="1433382370738" x="61.92553" y="305.7788" width="126.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc40e990c0a80c610191f5df36ddb9b2</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="66" label="is_a" layerID="1" created="1433382437823"
        x="187.06247" y="312.48416" width="177.96832" height="55.03589"
        strokeWidth="1.0" autoSized="false" controlCount="2"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc40e991c0a80c610191f5dfa3a027d6</URIString>
        <point1 x="187.56247" y="313.98056"/>
        <point2 x="364.5308" y="367.02005"/>
        <ID1 xsi:type="node">63</ID1>
        <ctrlPoint0 x="301.28674" y="307.99225" xsi:type="point"/>
        <ctrlPoint1 x="356.09595" y="329.07956" xsi:type="point"/>
    </child>
    <child ID="75" label="is_a" layerID="1" created="1433382954228"
        x="16.609653" y="273.8016" width="98.553024" height="29.366241"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc4a5783c0a80c610191f5df2f213378</URIString>
        <point1 x="18.109653" y="294.63504"/>
        <point2 x="113.662674" y="301.66785"/>
        <ID1 xsi:type="node">76</ID1>
        <ctrlPoint0 x="95.75846" y="268.45175" xsi:type="point"/>
    </child>
    <child ID="91" label="ncro: mammalian_miRNA" layerID="1"
        created="1433548846994" x="296.07248" y="369.75262"
        width="154.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c63bb3b9c0a80c4a00148dd84aae1949</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="100" label="is_a" layerID="1" created="1433549793675"
        x="408.64722" y="335.10516" width="313.27325" height="35.14746"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c63bb3bdc0a80c4a00148dd88ff8e05c</URIString>
        <point1 x="409.14722" y="369.75262"/>
        <point2 x="721.4205" y="335.60516"/>
        <ID1 xsi:type="node">91</ID1>
        <ID2 xsi:type="node">30</ID2>
        <ctrlPoint0 x="523.0944" y="333.42828" xsi:type="point"/>
    </child>
    <child ID="104" label="GO: transcription" layerID="1"
        created="1433550288241" x="868.6929" y="146.95157" width="103.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#B5B995</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c644475ac0a80c4a00148dd894a8124b</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="107" label="PRO: protein" layerID="1"
        created="1433550323639" x="107.582016" y="75.28398" width="81.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#83CEFF</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/c644475bc0a80c4a00148dd8352bab0c</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="131" label="RO: participates_in" layerID="1"
        created="1433750311138" x="925.412" y="169.45154" width="91.0"
        height="157.2284" strokeWidth="1.0" autoSized="false"
        controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2325cc2c0a80c4a00cfd279cd814c0c</URIString>
        <point1 x="1001.36115" y="326.17993"/>
        <point2 x="928.23126" y="169.95154"/>
        <ID1 xsi:type="node">46</ID1>
        <ID2 xsi:type="node">104</ID2>
        <ctrlPoint0 x="977.0277" y="239.76085" xsi:type="point"/>
    </child>
    <child ID="135" label="CLO: cell line" layerID="1"
        created="1433751135911" x="942.3611" y="-74.40573" width="84.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#BDE5F2</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d23e9359c0a80c4a00cfd2791b229471</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="140" label="BFO2: occurs_in" layerID="1"
        created="1433751458806" x="941.9882" y="-51.90576"
        width="100.086" height="199.35733" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d23f688cc0a80c4a00cfd2790718175d</URIString>
        <point1 x="942.4882" y="146.95157"/>
        <point2 x="988.1863" y="-51.40576"/>
        <ID1 xsi:type="node">104</ID1>
        <ID2 xsi:type="node">135</ID2>
        <ctrlPoint0 x="1037.8112" y="97.78384" xsi:type="point"/>
    </child>
    <child ID="143" label="ncro: is_model_of_disease" layerID="1"
        created="1433751632123" x="570.0147" y="-157.4559"
        width="399.43005" height="83.55014" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#7C18C9</strokeColor>
        <textColor>#7C18C9</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d24213edc0a80c4a00cfd279c290fa2e</URIString>
        <point1 x="968.94476" y="-74.40576"/>
        <point2 x="570.5147" y="-146.60599"/>
        <ID1 xsi:type="node">135</ID1>
        <ID2 xsi:type="node">126</ID2>
        <ctrlPoint0 x="819.0278" y="-186.23915" xsi:type="point"/>
    </child>
    <child ID="144" label="ncro: transcription_of_miRNA" layerID="1"
        created="1433751984266" x="447.6946" y="16.260838" width="172.0"
        height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d24bfccfc0a80c4a00cfd27909ff6308</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="145" label="is_a" layerID="1" created="1433752005899"
        x="559.8956" y="-2.1331916" width="350.01752" height="149.5848"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d24bfccfc0a80c4a00cfd27919c338fa</URIString>
        <point1 x="560.3956" y="16.260849"/>
        <point2 x="909.41315" y="146.9516"/>
        <ID1 xsi:type="node">144</ID1>
        <ID2 xsi:type="node">104</ID2>
        <ctrlPoint0 x="721.76105" y="-53.238384" xsi:type="point"/>
    </child>
    <child ID="147" label="ncro: miRNA_transcription_initiation"
        layerID="1" created="1433752577010" x="354.361" y="103.59421"
        width="208.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25561dac0a80c4a00cfd27933737d8f</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="148" label="ncro: protein_miRNA_promoter_binding"
        layerID="1" created="1433752594051" x="451.6944" y="193.59422"
        width="230.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25561dac0a80c4a00cfd279ddf26489</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="149" label="ncro: promoter_of_miRNA" layerID="1"
        created="1433752659308" x="645.02765" y="102.26087"
        width="154.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FDE888</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25561dac0a80c4a00cfd279322bbbff</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="151" label="RO: participate_in" layerID="1"
        created="1433753378695" x="184.07683" y="95.58069"
        width="268.86523" height="104.81018" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d25de720c0a80c4a00cfd27921d628da</URIString>
        <point1 x="184.57683" y="96.08069"/>
        <point2 x="452.44208" y="199.89087"/>
        <ID1 xsi:type="node">107</ID1>
        <ID2 xsi:type="node">148</ID2>
        <ctrlPoint0 x="434.23438" y="159.67857" xsi:type="point"/>
        <ctrlPoint1 x="349.56738" y="195.20567" xsi:type="point"/>
    </child>
    <child ID="152" label="ncro: regulate_miRNA_transcription"
        layerID="1" created="1433753551626" x="151.07207" y="-43.547184"
        width="380.63684" height="119.33119" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#AF55F4</strokeColor>
        <textColor>#AF55F4</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2615c64c0a80c4a00cfd27950d3405c</URIString>
        <point1 x="151.57207" y="75.284"/>
        <point2 x="531.2089" y="16.260742"/>
        <ID1 xsi:type="node">107</ID1>
        <ID2 xsi:type="node">144</ID2>
        <ctrlPoint0 x="201.49495" y="-89.215996" xsi:type="point"/>
        <ctrlPoint1 x="519.31714" y="-38.758076" xsi:type="point"/>
    </child>
    <child ID="154" label="RO: precedes" layerID="1"
        created="1433753916773" x="447.03973" y="126.09424"
        width="92.138306" height="67.999985" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d26514fcc0a80c4a00cfd279d6be07aa</URIString>
        <point1 x="538.67804" y="193.59422"/>
        <point2 x="459.48517" y="126.59424"/>
        <ID1 xsi:type="node">148</ID1>
        <ID2 xsi:type="node">147</ID2>
        <ctrlPoint0 x="462.99786" y="162.52948" xsi:type="point"/>
    </child>
    <child ID="155" label="RO: precedes" layerID="1"
        created="1433753966436" x="424.11267" y="38.760742"
        width="107.6203" height="64.18831" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2659098c0a80c4a00cfd279740792b6</URIString>
        <point1 x="424.61267" y="102.44905"/>
        <point2 x="531.233" y="39.260742"/>
        <ID2 xsi:type="node">144</ID2>
        <ctrlPoint0 x="455.01517" y="25.4114" xsi:type="point"/>
        <ctrlPoint1 x="517.1143" y="105.21695" xsi:type="point"/>
    </child>
    <child ID="156" label="RO: participates_in" layerID="1"
        created="1433754046535" x="619.1946" y="27.309229"
        width="150.83728" height="75.45163" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd279f89efdba</URIString>
        <point1 x="726.40967" y="102.260864"/>
        <point2 x="619.6946" y="27.809229"/>
        <ID1 xsi:type="node">149</ID1>
        <ID2 xsi:type="node">144</ID2>
        <ctrlPoint0 x="742.086" y="61.121044" xsi:type="point"/>
        <ctrlPoint1 x="741.29755" y="27.877655" xsi:type="point"/>
    </child>
    <child ID="158" label="GO: RNA_metabolic_process" layerID="1"
        created="1433754139438" x="229.19443" y="24.594929"
        width="173.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#B5B995</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd2792b44e1da</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="159" label="is_a" layerID="1" created="1433754164170"
        x="316.11972" y="47.09497" width="67.1673" height="56.999237"
        strokeWidth="1.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#404040</strokeColor>
        <textColor>#404040</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd2793c92f5c8</URIString>
        <point1 x="382.78702" y="103.59421"/>
        <point2 x="316.61972" y="47.59497"/>
        <ID1 xsi:type="node">147</ID1>
        <ID2 xsi:type="node">158</ID2>
        <ctrlPoint0 x="320.36108" y="94.09493" xsi:type="point"/>
    </child>
    <child ID="164" label="Hox-D10 gene" layerID="1"
        created="1433756972501" x="647.0742" y="616.9813" width="112.0"
        height="27.0" strokeWidth="3.0" autoSized="true" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d298ed21c0a80c4a008a5fd2a49200fc</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="165" label="is_a" layerID="1" created="1433757175291"
        x="690.00415" y="546.8197" width="30.0" height="71.66174"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d298ed21c0a80c4a008a5fd25978d7e5</URIString>
        <point1 x="704.0358" y="616.98145"/>
        <point2 x="701.8325" y="548.3197"/>
        <ID1 xsi:type="node">164</ID1>
        <ctrlPoint0 x="707.07416" y="574.3253" xsi:type="point"/>
    </child>
    <child ID="76" label="hsa-miR-10b" layerID="1"
        created="1433382974053" x="-72.42582" y="294.63504"
        width="101.0" height="27.0" strokeWidth="3.0" autoSized="true" xsi:type="node">
        <fillColor>#00000000</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/bc4a5782c0a80c610191f5dfd082b240</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="167" label="MDA-MB-231 cell" layerID="1"
        created="1433758038049" x="1048.5742" y="16.457031"
        width="143.16748" height="27.0" strokeWidth="3.0"
        autoSized="false" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a4aec2c0a80c4a008a5fd2c0922196</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="168" label="is_a" layerID="1" created="1433758038049"
        x="1028.1476" y="-63.485527" width="92.66223" height="81.44256"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a4aec2c0a80c4a008a5fd293399dec</URIString>
        <point1 x="1119.3098" y="16.457031"/>
        <point2 x="1029.6476" y="-61.985527"/>
        <ID1 xsi:type="node">167</ID1>
        <ctrlPoint0 x="1116.0742" y="-35.042953" xsi:type="point"/>
    </child>
    <child ID="170" label="is_a" layerID="1" created="1433758161095"
        x="300.0272" y="-170.97076" width="180.897" height="24.104965"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a7b360c0a80c4a008a5fd2716b7de7</URIString>
        <point1 x="301.5272" y="-148.3658"/>
        <point2 x="479.4242" y="-148.41908"/>
        <ID1 xsi:type="node">171</ID1>
        <ctrlPoint0 x="368.26758" y="-174.54909" xsi:type="point"/>
    </child>
    <child ID="171" label="metastatic breast cancer" layerID="1"
        created="1433758161095" x="168.14886" y="-148.3658"
        width="197.93445" height="27.0" strokeWidth="3.0"
        autoSized="false" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a7b360c0a80c4a008a5fd21a1e5ad7</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="126" label="DOID: disease" layerID="1"
        created="1433715751553" x="480.40753" y="-151.07172"
        width="92.0" height="23.0" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#C1F780</fillColor>
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d01f8535c0a80c4a019878a268cc4552</URIString>
        <shape arcwidth="20.0" archeight="20.0" xsi:type="roundRect"/>
    </child>
    <child ID="172" label="TWIST1" layerID="1" created="1433758301418"
        x="62.017296" y="176.19014" width="67.0" height="27.0"
        strokeWidth="3.0" autoSized="true" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a7b36ac0a80c4a008a5fd22ea8013f</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="173" label="is_a" layerID="1" created="1433758301418"
        x="108.47801" y="98.84091" width="47.13829" height="78.84923"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a7b36ac0a80c4a008a5fd221f3d997</URIString>
        <point1 x="109.97801" y="176.19014"/>
        <point2 x="146.33878" y="100.34091"/>
        <ID1 xsi:type="node">172</ID1>
        <ctrlPoint0 x="153.07422" y="135.95705" xsi:type="point"/>
    </child>
    <child ID="176" label="hsa-miR-10b transcription" layerID="1"
        created="1433758386518" x="515.96155" y="-79.33702"
        width="201.11267" height="27.0" strokeWidth="3.0"
        autoSized="false" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a92f12c0a80c4a008a5fd268cc9a6f</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="177" label="is_a" layerID="1" created="1433758386518"
        x="544.33044" y="-53.837036" width="63.661255" height="67.13625"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2a92f12c0a80c4a008a5fd2fd202014</URIString>
        <point1 x="606.4917" y="-52.337036"/>
        <point2 x="545.83044" y="11.799219"/>
        <ID1 xsi:type="node">176</ID1>
        <ctrlPoint0 x="586.5742" y="-25.518661" xsi:type="point"/>
    </child>
    <child ID="178" label="TWIST1 binding to hsa-miR-10b E-box 1"
        layerID="1" created="1433758788450" x="343.41455" y="288.81937"
        width="304.15173" height="27.0" strokeWidth="3.0"
        autoSized="false" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b29d39c0a80c4a00ffd2a35f2f5fea</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="179" label="is_a" layerID="1" created="1433758788450"
        x="506.32248" y="217.95058" width="39.999176" height="72.36876"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b29d39c0a80c4a00ffd2a3a47b0207</URIString>
        <point1 x="507.82248" y="288.81934"/>
        <point2 x="540.5601" y="219.45058"/>
        <ID1 xsi:type="node">178</ID1>
        <ctrlPoint0 x="538.45197" y="255.28888" xsi:type="point"/>
    </child>
    <child ID="180" label="hsa-miR-10b transcription initiation"
        layerID="1" created="1433758918848" x="142.32388" y="215.44391"
        width="267.10043" height="27.0" strokeWidth="3.0"
        autoSized="false" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b29d43c0a80c4a00ffd2a3c9e87d87</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="181" label="is_a" layerID="1" created="1433758918848"
        x="225.6865" y="112.54379" width="126.64142" height="104.400116"
        strokeWidth="3.0" autoSized="false" controlCount="2"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b29d43c0a80c4a00ffd2a31800540c</URIString>
        <point1 x="256.69644" y="215.44391"/>
        <point2 x="350.8279" y="114.04379"/>
        <ID1 xsi:type="node">180</ID1>
        <ctrlPoint0 x="224.01529" y="192.43819" xsi:type="point"/>
        <ctrlPoint1 x="215.30722" y="113.83115" xsi:type="point"/>
    </child>
    <child ID="182" label="hsa-miR-10b E-box 1" layerID="1"
        created="1433759075936" x="737.46814" y="-73.45655"
        width="159.0" height="27.0" strokeWidth="3.0" autoSized="true" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b530b7c0a80c4a00ffd2a3c1979e6e</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="183" label="is_a" layerID="1" created="1433759075936"
        x="793.09045" y="-47.956543" width="59.33606" height="149.79745"
        strokeWidth="3.0" autoSized="false" controlCount="2"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b530b7c0a80c4a00ffd2a3e343386a</URIString>
        <point1 x="815.1924" y="-46.456543"/>
        <point2 x="794.59045" y="100.34091"/>
        <ID1 xsi:type="node">182</ID1>
        <ctrlPoint0 x="807.07544" y="15.251873" xsi:type="point"/>
        <ctrlPoint1 x="889.4677" y="26.82088" xsi:type="point"/>
    </child>
    <child ID="184" label="Hox-D10 transcription" layerID="1"
        created="1433759166587" x="1049.5516" y="204.00133"
        width="169.34668" height="27.0" strokeWidth="3.0"
        autoSized="false" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b530b7c0a80c4a00ffd2a365233fc3</URIString>
        <shape xsi:type="rectangle"/>
    </child>
    <child ID="185" label="is_a" layerID="1" created="1433759166587"
        x="972.282" y="149.77092" width="160.849" height="55.730545"
        strokeWidth="3.0" autoSized="false" controlCount="1"
        arrowState="2" xsi:type="link">
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-plain-16</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2b530c1c0a80c4a00ffd2a32d33a202</URIString>
        <point1 x="1131.631" y="204.00146"/>
        <point2 x="973.782" y="155.15244"/>
        <ID1 xsi:type="node">184</ID1>
        <ctrlPoint0 x="1119.1342" y="138.96489" xsi:type="point"/>
    </child>
    <child ID="186"
        label="Ma L, Teruya-Feldstein J, Weinberg RA. Tumour invasion and metastasis&#xa;initiated by microRNA-10b in breast cancer. Nature. 2007 Oct 11;449(7163):682-8. &#xa;Epub 2007 Sep 26. Erratum in: Nature. 2008 Sep 11;455(7210):256.&#xa;PubMed PMID: 17898713."
        layerID="1" created="1433761334951" x="-46.42578" y="530.62915"
        width="525.8" height="116.8" strokeWidth="1.0" autoSized="true" xsi:type="node">
        <fillColor>#FFFFFF</fillColor>
        <strokeColor>#DD7B11</strokeColor>
        <textColor>#DD7B11</textColor>
        <font>Arial-bold-12</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2d81d24c0a80c4a01dc840ba544c775</URIString>
        <shape xsi:type="octagon"/>
    </child>
    <child ID="187" label="RO: participates_in" layerID="1"
        created="1433762634398" x="555.23413" y="123.26686"
        width="95.35077" height="70.82738" strokeWidth="1.0"
        autoSized="false" controlCount="1" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2ead2dac0a80c4a0073f7fd2991206b</URIString>
        <point1 x="650.0849" y="123.76686"/>
        <point2 x="570.54095" y="193.59424"/>
        <ID1 xsi:type="node">149</ID1>
        <ID2 xsi:type="node">148</ID2>
        <ctrlPoint0 x="591.1554" y="131.96295" xsi:type="point"/>
    </child>
    <child ID="157" label="RO: participates_in" layerID="1"
        created="1433754091599" x="540.0601" y="57.505966"
        width="163.32782" height="45.443085" strokeWidth="1.0"
        autoSized="false" controlCount="2" arrowState="2" xsi:type="link">
        <strokeColor>#000000</strokeColor>
        <textColor>#000000</textColor>
        <font>Arial-plain-11</font>
        <URIString>http://vue.tufts.edu/rdf/resource/d2693664c0a80c4a00cfd2795daa83d2</URIString>
        <point1 x="702.88794" y="102.260864"/>
        <point2 x="540.5601" y="102.44905"/>
        <ID1 xsi:type="node">149</ID1>
        <ctrlPoint0 x="629.799" y="58.345726" xsi:type="point"/>
        <ctrlPoint1 x="540.72327" y="45.43354" xsi:type="point"/>
    </child>
    <layer ID="1" label="Layer 1" created="1431734769931" x="0.0"
        y="0.0" width="1.4E-45" height="1.4E-45" strokeWidth="0.0" autoSized="false">
        <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b60a147a5701bb3ae4dfffb6d3</URIString>
    </layer>
    <userZoom>1.0</userZoom>
    <userOrigin x="-87.92582" y="-295.13528"/>
    <presentationBackground>#202020</presentationBackground>
    <PathwayList currentPathway="0" revealerIndex="-1">
        <pathway ID="0" label="Untitled Pathway" created="1431734769931"
            x="0.0" y="0.0" width="1.4E-45" height="1.4E-45"
            strokeWidth="0.0" autoSized="false" currentIndex="0" open="true">
            <strokeColor>#B3CC33CC</strokeColor>
            <textColor>#000000</textColor>
            <font>SansSerif-plain-14</font>
            <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b70a147a5701bb3ae4c573cfec</URIString>
            <masterSlide ID="2" created="1431734769956" x="0.0" y="0.0"
                width="800.0" height="600.0" locked="true"
                strokeWidth="0.0" autoSized="false">
                <fillColor>#000000</fillColor>
                <strokeColor>#404040</strokeColor>
                <textColor>#000000</textColor>
                <font>SansSerif-plain-14</font>
                <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b80a147a5701bb3ae4de1adb6c</URIString>
                <titleStyle ID="3" label="Header"
                    created="1431734769971" x="335.5" y="171.0"
                    width="129.0" height="58.0" strokeWidth="0.0"
                    autoSized="true" isStyle="true" xsi:type="node">
                    <strokeColor>#404040</strokeColor>
                    <textColor>#FFFFFF</textColor>
                    <font>Gill Sans-plain-36</font>
                    <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b80a147a5701bb3ae4d4c9c2ee</URIString>
                    <shape xsi:type="rectangle"/>
                </titleStyle>
                <textStyle ID="4" label="Slide Text"
                    created="1431734769972" x="346.5" y="280.5"
                    width="107.0" height="39.0" strokeWidth="0.0"
                    autoSized="true" isStyle="true" xsi:type="node">
                    <strokeColor>#404040</strokeColor>
                    <textColor>#FFFFFF</textColor>
                    <font>Gill Sans-plain-22</font>
                    <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b90a147a5701bb3ae4b6db12b2</URIString>
                    <shape xsi:type="rectangle"/>
                </textStyle>
                <linkStyle ID="5" label="Links" created="1431734769973"
                    x="373.5" y="383.0" width="53.0" height="34.0"
                    strokeWidth="0.0" autoSized="true" isStyle="true" xsi:type="node">
                    <strokeColor>#404040</strokeColor>
                    <textColor>#B3BFE3</textColor>
                    <font>Gill Sans-plain-18</font>
                    <URIString>http://vue.tufts.edu/rdf/resource/5a6f92b90a147a5701bb3ae44e67e965</URIString>
                    <shape xsi:type="rectangle"/>
                </linkStyle>
            </masterSlide>
        </pathway>
    </PathwayList>
    <date>2015-05-15</date>
    <modelVersion>6</modelVersion>
    <saveLocation>C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design</saveLocation>
    <saveFile>C:\Users\Jingshan Huang\Documents\GitHub\OmniSearch\ncRO-Ontology-Design\ncRO-UseCase-1-step-1.vue</saveFile>
</LW-MAP>
