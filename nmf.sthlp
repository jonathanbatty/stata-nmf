{smcl}
{* *! version 1.0 15 Sep 2022}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install command2" "ssc install command2"}{...}
{vieweralsosee "Help command2 (if installed)" "help command2"}{...}
{viewerjumpto "Syntax" "nmf##syntax"}{...}
{viewerjumpto "Description" "nmf##description"}{...}
{viewerjumpto "Options" "nmf##options"}{...}
{viewerjumpto "Remarks" "nmf##remarks"}{...}
{viewerjumpto "Examples" "nmf##examples"}{...}
{title:Title}
{phang}
{bf:nmf} {hline 2} <Insert title>

{marker syntax}{...}
{title:Syntax}
{p 8 17 2}
{cmdab:nmf}
varlist(numeric)
[{cmd:,}
{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}

{syntab:Required }
{synopt:{opt k(#)}} {synopt:{opt iter(#)}} 
{syntab:Optional}
{synopt:{opt initial(string)}}  

{synopt:{opt method(string)}}  

{synopt:{opt loss(string)}}  

{synopt:{opt stop(numlist max  =  1)}}  

{synopt:{opt nograph}}  

{synopt:{opt noframes}}  

{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}
{pstd}

{marker options}{...}
{title:Options}
{dlgtab:Main}

{phang}
{opt k(#)}  

{phang}
{opt iter(#)}  

{phang}
{opt initial(string)}  

{phang}
{opt method(string)}  

{phang}
{opt loss(string)}  

{phang}
{opt stop(numlist max  =  1)}  

{phang}
{opt nograph}  

{phang}
{opt noframes}  



{marker examples}{...}
{title:Examples}

{title:Stored results}

{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:r(W)}}  {p_end}
{synopt:{cmd:r(H)}}  {p_end}
{synopt:{cmd:r(norms)}}  {p_end}


{title:Author}
{p}



