// .vuepress/config.js
module.exports = {
	locales: {
		'/': {
			lang: 'en-US',
			title: 'ReacNetGenerator',
			description: 'An automatic generator of reaction network for reactive molecular dynamics simulation'
		},
		'/zh/': {
			lang: 'zh-CN',
			title: 'ReacNetGenerator',
			description: '反应动力学模拟的反应网络自动生成器'
		}
	},
	head: [
		['link', { rel: 'icon', href: `/reacnetgen.svg` }],
		['meta', { name: 'apple-mobile-web-app-capable', content: 'yes' }],
		['link', { rel: 'apple-touch-icon', href: `/reacnetgen.svg` }],
		['meta', { name: 'msapplication-TileImage', content: '/reacnetgen.svg' }]
	],
	themeConfig: {
		repo: "tongzhugroup/reacnetgenerator",
		editLinks: true,
		docsDir: 'docs',
		smoothScroll: true,
		locales: {
			'/': {
				selectText: 'Languages',
				label: 'English',
				nav: [
					{ text: 'Home', link: '/' },
					{ text: 'Guide', link: '/guide/' },
					{ text: 'Tutorial', link: '/tutorial/' },
					{
						text: 'See more',
						ariaLabel: 'See more',
						items: [
							{ text: 'Report', link: '/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json&', target: '_blank' },
							{ text: 'Article', link: 'https://doi.org/10.1039/C9CP05091D' },
							{ text: 'Group', link: 'https://computchem.cn' },
						]
					}
				],
				sidebar: {
					'/guide/': getGuideSidebar("Guide"),
					'/tutorial/': getTutorialSidebar("Tutorial"),
				},
			},
			'/zh/': {
				selectText: '语言',
				label: '中文',
				nav: [
					{ text: '主页', link: '/zh/' },
					{ text: '指南', link: '/zh/guide/' },
					{ text: '教程', link: '/tutorial/' },
					{
						text: '了解更多',
						ariaLabel: '了解更多',
						items: [
							{ text: '分析结果', link: '/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json&', target: '_blank' },
							{ text: '论文', link: 'https://doi.org/10.1039/C9CP05091D' },
							{ text: '课题组', link: 'https://computchem.cn' },
						]
					}
				],
				sidebar: {
					'/zh/guide/': getGuideSidebar("指南"),
					'/zh/tutorial/': getTutorialSidebar("教程"),
				},
			},
		}
	},
}


function getGuideSidebar(title) {
	return [
		{
			title: title,
			collapsable: true,
			children: [
				"",
				"build",
			]
		}
	]
}


function getTutorialSidebar(title) {
	return [
		{
			title: title,
			collapsable: true,
			children: [
				"",
				"install",
			]
		}
	]
}
