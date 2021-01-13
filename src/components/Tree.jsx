export const RuleTree = props => {
  return (
    <ul>
      {props.rules.map(rule => {
        if (!rule.children)
          return (<li key={rule.uuid}>{rule.name}</li>)
        return (<li key={rule.uuid}>{rule.name}<RuleTree rules={rule.children}/></li>)
      })}
    </ul>
  )
}
