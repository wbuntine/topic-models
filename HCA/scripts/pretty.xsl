<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
  <xsl:output method="text" encoding="UTF-8"/>
  <xsl:template match="text()"/>
  <xsl:template match="//Documents">
    <xsl:text>Documents </xsl:text><xsl:value-of select="normalize-space(@count)"/>
    <xsl:text>
</xsl:text>
  </xsl:template>
  <xsl:template match="//Topic">
    <xsl:text>Topic </xsl:text><xsl:value-of select="@id"/><xsl:text> (</xsl:text><xsl:value-of select="@prob"/><xsl:text>):</xsl:text>
    <xsl:for-each select='Word[position()&lt;8]'>
    <xsl:text> </xsl:text><xsl:value-of select="."/><xsl:text>,</xsl:text>
    </xsl:for-each>
    <xsl:text>
</xsl:text>
  </xsl:template>
</xsl:stylesheet>



